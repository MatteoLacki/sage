use super::Batch;
use anyhow::{ensure, Result};
use sage_core::intensity_codec::IntensityEncoding;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub struct Writer {
    root: PathBuf,
    intensity: BufWriter<File>,
    annotation: BufWriter<File>,
    starts: BufWriter<File>,
    ends: BufWriter<File>,
    values: Vec<u8>,
    ids: Vec<u8>,
    encoding: IntensityEncoding,
    pub fragments: i64,
    pub slots: usize,
}

fn file(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        1024 * 1024,
        OpenOptions::new().write(true).create_new(true).open(path)?,
    ))
}

impl Writer {
    pub fn new(root: &Path, encoding: IntensityEncoding) -> Result<Self> {
        std::fs::create_dir(root)?;
        let arrays = root.join("arrays.mmappet");
        let index = root.join("index.mmappet");
        std::fs::create_dir(&arrays)?;
        std::fs::create_dir(&index)?;
        std::fs::write(
            arrays.join("schema.txt"),
            match encoding {
                IntensityEncoding::Fp16 => "float16 predicted_intensity\nuint8 annotation_id\n",
                IntensityEncoding::Sqrt8Basepeak => {
                    "uint8 predicted_intensity\nuint8 annotation_id\n"
                }
            },
        )?;
        if encoding == IntensityEncoding::Sqrt8Basepeak {
            std::fs::write(
                arrays.join("meta.json"),
                r#"{"version":2,"intensity_encoding":"sqrt8_basepeak_v1","normalization":"per_spectrum_base_peak","rounding":"nearest_even"}"#,
            )?;
        }
        std::fs::write(index.join("schema.txt"), "int64 start\nint64 end\n")?;
        Ok(Self {
            root: root.to_owned(),
            intensity: file(&arrays.join("0.bin"))?,
            annotation: file(&arrays.join("1.bin"))?,
            starts: file(&index.join("0.bin"))?,
            ends: file(&index.join("1.bin"))?,
            values: Vec::new(),
            ids: Vec::new(),
            encoding,
            fragments: 0,
            slots: 0,
        })
    }

    pub fn append(&mut self, batch: &Batch) -> Result<()> {
        ensure!(
            batch.first_slot == self.slots,
            "out-of-order prediction batch"
        );
        self.values.clear();
        self.ids.clear();
        let mut valid = 0;
        for slot in 0..batch.slots {
            let Some(_) = batch.inputs.get(valid).filter(|input| input.slot == slot) else {
                self.starts.write_all(&(-1_i64).to_le_bytes())?;
                self.ends.write_all(&(-1_i64).to_le_bytes())?;
                continue;
            };
            let start = self.fragments;
            let (values, ids) = batch.output.row(valid)?;
            let expected_width = if self.encoding == IntensityEncoding::Fp16 {
                2
            } else {
                1
            };
            ensure!(
                batch.output.width() == expected_width,
                "prediction encoding mismatch"
            );
            self.values.extend_from_slice(values);
            self.ids.extend_from_slice(ids);
            self.fragments = self
                .fragments
                .checked_add(ids.len() as i64)
                .ok_or_else(|| anyhow::anyhow!("fragment offset overflow"))?;
            self.starts.write_all(&start.to_le_bytes())?;
            self.ends.write_all(&self.fragments.to_le_bytes())?;
            valid += 1;
        }
        ensure!(
            valid == batch.inputs.len(),
            "unconsumed predictions in batch"
        );
        self.intensity.write_all(&self.values)?;
        self.annotation.write_all(&self.ids)?;
        self.slots += batch.slots;
        Ok(())
    }

    pub fn finish(mut self) -> Result<(usize, i64)> {
        for writer in [
            &mut self.intensity,
            &mut self.annotation,
            &mut self.starts,
            &mut self.ends,
        ] {
            writer.flush()?;
            writer.get_ref().sync_all()?;
        }
        for path in ["arrays.mmappet/schema.txt", "index.mmappet/schema.txt"] {
            File::open(self.root.join(path))?.sync_all()?;
        }
        if self.encoding == IntensityEncoding::Sqrt8Basepeak {
            File::open(self.root.join("arrays.mmappet/meta.json"))?.sync_all()?;
        }
        Ok((self.slots, self.fragments))
    }
}
