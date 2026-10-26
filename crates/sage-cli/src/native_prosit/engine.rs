//! A dynamically loaded C ABI keeps TensorRT/CUDA out of ordinary Sage builds.
use anyhow::{bail, ensure, Result};
use sage_core::intensity_codec::IntensityEncoding;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

#[link(name = "dl")]
extern "C" {
    fn dlopen(filename: *const c_char, flags: c_int) -> *mut c_void;
    fn dlsym(handle: *mut c_void, symbol: *const c_char) -> *mut c_void;
    fn dlerror() -> *const c_char;
    fn dlclose(handle: *mut c_void) -> c_int;
}

/// ABI v2: one count byte, K intensity values, then K annotation IDs per row.
/// Unused entries are padding, never written to the library. Count 255 is an error.
pub struct Predictions {
    pub(crate) data: Vec<u8>,
    k: usize,
    encoding: IntensityEncoding,
}

impl Predictions {
    pub fn new(top_k: usize, encoding: IntensityEncoding) -> Self {
        Self {
            data: Vec::new(),
            k: if top_k == 0 { 174 } else { top_k.min(174) },
            encoding,
        }
    }
    pub(crate) fn width(&self) -> usize {
        if self.encoding == IntensityEncoding::Fp16 {
            2
        } else {
            1
        }
    }
    pub(crate) fn stride(&self) -> usize {
        1 + self.k * (self.width() + 1)
    }
    pub(crate) fn resize(&mut self, rows: usize) {
        self.data.resize(rows * self.stride(), 0);
    }
    pub(crate) fn row(&self, row: usize) -> Result<(&[u8], &[u8])> {
        let stride = self.stride();
        let record = self
            .data
            .get(row * stride..(row + 1) * stride)
            .ok_or_else(|| anyhow::anyhow!("missing packed prediction row {row}"))?;
        let count = record[0] as usize;
        ensure!(
            count != 255,
            "nonfinite model output, FP16 overflow, or invalid GPU input at prediction {row}"
        );
        ensure!(count <= self.k, "invalid packed fragment count {count}");
        let ids = 1 + self.k * self.width();
        Ok((
            &record[1..1 + count * self.width()],
            &record[ids..ids + count],
        ))
    }
}

pub trait Engine {
    fn predict(
        &mut self,
        tokens: &[i32],
        charges: &[f32],
        energies: &[f32],
        output: &mut Predictions,
    ) -> Result<()>;
}

type Create = unsafe extern "C" fn(*const c_char, c_int, *mut c_char, usize) -> *mut c_void;
type Infer = unsafe extern "C" fn(
    *mut c_void,
    c_int,
    *const i32,
    *const f32,
    *const f32,
    c_int,
    c_int,
    *mut u8,
    *mut c_char,
    usize,
) -> c_int;
type Destroy = unsafe extern "C" fn(*mut c_void);

struct Library(*mut c_void);
impl Drop for Library {
    fn drop(&mut self) {
        unsafe {
            dlclose(self.0);
        }
    }
}
impl Library {
    fn symbol(&self, name: &str) -> Result<*mut c_void> {
        let name = CString::new(name)?;
        let pointer = unsafe { dlsym(self.0, name.as_ptr()) };
        ensure!(!pointer.is_null(), "missing native symbol {:?}", name);
        Ok(pointer)
    }
}

pub struct TensorRt {
    // Destroy the model while the library is still loaded.
    _library: Library,
    model: *mut c_void,
    infer: Infer,
    destroy: Destroy,
    capacity: usize,
    output_kind: u32,
}

impl TensorRt {
    pub fn open(library: &Path, plan: &Path, capacity: usize) -> Result<Self> {
        ensure!((1..=1000).contains(&capacity), "batch size must be 1..1000");
        let filename = CString::new(library.as_os_str().as_bytes())?;
        let pointer = unsafe { dlopen(filename.as_ptr(), 2) }; // RTLD_NOW | RTLD_LOCAL
        if pointer.is_null() {
            let error = unsafe { dlerror() };
            let message = if error.is_null() {
                "unknown dlopen error".into()
            } else {
                unsafe { CStr::from_ptr(error) }
                    .to_string_lossy()
                    .into_owned()
            };
            bail!("load {}: {message}", library.display());
        }
        let library = Library(pointer);
        // SAFETY: bridge ABI version and symbol signatures are defined in prosit_trt.cpp.
        let version: unsafe extern "C" fn() -> u32 =
            unsafe { std::mem::transmute(library.symbol("sage_prosit_abi_version")?) };
        ensure!(
            unsafe { version() } == 3,
            "unsupported native bridge ABI; rebuild the CUDA bridge (requires ABI 3)"
        );
        let create: Create = unsafe { std::mem::transmute(library.symbol("sage_prosit_create")?) };
        let infer: Infer = unsafe { std::mem::transmute(library.symbol("sage_prosit_infer")?) };
        let destroy: Destroy =
            unsafe { std::mem::transmute(library.symbol("sage_prosit_destroy")?) };
        let output_kind: unsafe extern "C" fn(*mut c_void) -> u32 =
            unsafe { std::mem::transmute(library.symbol("sage_prosit_output_kind")?) };
        let filename = CString::new(plan.as_os_str().as_bytes())?;
        let mut error = [0 as c_char; 4096];
        let model = unsafe {
            create(
                filename.as_ptr(),
                capacity as c_int,
                error.as_mut_ptr(),
                error.len(),
            )
        };
        ensure!(
            !model.is_null(),
            "TensorRT initialization: {}",
            unsafe { CStr::from_ptr(error.as_ptr()) }.to_string_lossy()
        );
        let output_kind = unsafe { output_kind(model) };
        Ok(Self {
            _library: library,
            model,
            infer,
            destroy,
            capacity,
            output_kind,
        })
    }

    pub fn dense_output(&self) -> bool {
        self.output_kind != 0
    }

    pub fn model_output(&self) -> &'static str {
        match self.output_kind {
            1 => "dense_fp16_29x6",
            2 => "dense_fp32_29x6",
            _ => "legacy_post_activation_fp16_174",
        }
    }
}

impl Engine for TensorRt {
    fn predict(
        &mut self,
        tokens: &[i32],
        charges: &[f32],
        energies: &[f32],
        output: &mut Predictions,
    ) -> Result<()> {
        let n = energies.len();
        ensure!(
            n > 0
                && n <= self.capacity
                && tokens.len() == n * 30
                && charges.len() == n * 6
                && output.data.len() == n * output.stride(),
            "invalid native input/output buffer lengths"
        );
        let mut error = [0 as c_char; 4096];
        // All buffers remain alive until the bridge synchronizes its CUDA stream.
        let status = unsafe {
            (self.infer)(
                self.model,
                n as c_int,
                tokens.as_ptr(),
                charges.as_ptr(),
                energies.as_ptr(),
                output.k as c_int,
                if output.encoding == IntensityEncoding::Fp16 {
                    0
                } else {
                    1
                },
                output.data.as_mut_ptr(),
                error.as_mut_ptr(),
                error.len(),
            )
        };
        ensure!(
            status == 0,
            "TensorRT inference: {}",
            unsafe { CStr::from_ptr(error.as_ptr()) }.to_string_lossy()
        );
        Ok(())
    }
}

impl Drop for TensorRt {
    fn drop(&mut self) {
        unsafe {
            (self.destroy)(self.model);
        }
    }
}
