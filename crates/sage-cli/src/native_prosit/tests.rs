use super::*;
use sage_cloudpath::fragment_intensity_cache::{
    read_fragment_intensity_index, FragmentIntensityArrays,
};
use sage_core::peptide_view::PeptideView;
use sage_core::{database::Builder, fasta::Fasta, prosit::ModificationPolicy};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

fn directory() -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "sage_native_prosit_{}_{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

fn table() -> CompactPeptides {
    let fasta = Fasta::parse(">test\nMACDMKRACDMK\n".into(), "rev_", true);
    let mut builder: Builder = serde_json::from_value(serde_json::json!({
        "enzyme": {"cleave_at": "", "restrict": null, "min_len": 5, "max_len": 9},
        "peptide_min_mass": 0.0, "peptide_max_mass": 100000.0,
        "static_mods": {}, "variable_mods": {"[": [42.010565]}, "max_variable_mods": 1
    }))
    .unwrap();
    builder.fasta = Some("fixture".into());
    CompactPeptides::build(&builder.make_parameters(), &fasta).unwrap()
}

struct CheckedEngine {
    calls: usize,
    fail: bool,
}
impl engine::Engine for CheckedEngine {
    fn predict(
        &mut self,
        tokens: &[i32],
        charges: &[f32],
        energies: &[f32],
        output: &mut engine::Predictions,
    ) -> Result<()> {
        self.calls += 1;
        ensure!(!self.fail, "injected GPU failure");
        assert_eq!(tokens.len(), energies.len() * 30);
        assert_eq!(output.data.len(), energies.len() * output.stride());
        for input in tokens.chunks_exact(30) {
            let length = input.iter().position(|&t| t == 0).unwrap_or(30);
            assert!((5..=9).contains(&length));
            assert!(input[..length].iter().all(|&t| (1..=21).contains(&t)));
            assert!(input[length..].iter().all(|&t| t == 0));
        }
        for charge in charges.chunks_exact(6) {
            assert_eq!(charge.iter().sum::<f32>(), 1.0);
            assert!(charge[1] == 1.0 || charge[2] == 1.0);
        }
        assert!(energies.iter().all(|&e| e == 30.0_f32 * 0.01_f32));
        // Fake engine returns packed equal intensities; GPU correctness has a
        // separate CUDA fixture test. This exercises pipeline ordering and I/O.
        let width = output.width();
        let stride = output.stride();
        let k = (stride - 1) / (width + 1);
        for row in 0..energies.len() {
            let length = tokens[row * 30..(row + 1) * 30]
                .iter()
                .position(|&t| t == 0)
                .unwrap_or(30);
            let charge = charges[row * 6..(row + 1) * 6]
                .iter()
                .position(|&z| z == 1.0)
                .unwrap()
                + 1;
            let record = &mut output.data[row * stride..(row + 1) * stride];
            let mut count = 0;
            for id in 0..174 {
                let ordinal = id / 6 + 1;
                let z = id % 3;
                if ordinal >= length || z >= charge || count == k {
                    continue;
                }
                if width == 2 {
                    record[1 + count * 2..3 + count * 2].copy_from_slice(&0x3c00_u16.to_le_bytes());
                } else {
                    record[1 + count] = 255;
                }
                record[1 + k * width + count] = if id % 6 < 3 {
                    87 + (length - 1 - ordinal) * 3 + z
                } else {
                    (ordinal - 1) * 3 + z
                } as u8;
                count += 1;
            }
            record[0] = count as u8;
        }
        Ok(())
    }
}

#[test]
fn compact_pipeline_round_trips_skips_and_charge_slots() {
    let table = table();
    let tokenizer = Tokenizer::new(ModificationPolicy::Original, HashMap::new());
    let options = Options {
        row_start: 2,
        row_count: table.rows.len() - 3,
        min_charge: 2,
        max_charge: 3,
        collision_energy: 30.0,
        batch_size: 7,
        top_k: 0,
        intensity_encoding: IntensityEncoding::Fp16,
        skip_unsupported: true,
    };
    let path = directory();
    let mut engine = CheckedEngine {
        calls: 0,
        fail: false,
    };
    let metrics = run(&table, &tokenizer, &mut engine, &options, &path).unwrap();
    std::fs::write(
        path.join("index.mmappet/meta.json"),
        r#"{"dumped_peptides_sha256":"test"}"#,
    )
    .unwrap();
    let (fingerprint, index) = read_fragment_intensity_index(&path.join("index.mmappet")).unwrap();
    let arrays = FragmentIntensityArrays::open(&path.join("arrays.mmappet")).unwrap();
    assert_eq!(fingerprint, "test");
    assert_eq!(index.start.len(), options.row_count * 2);
    let mut cursor = 0;
    let mut skipped = 0;
    for slot in 0..index.start.len() {
        let row = options.row_start + slot / 2;
        let view = table.view(row);
        let mut tokens = [0; 30];
        if tokenizer.encode(&view, &mut tokens).is_err() {
            assert_eq!(index.start.as_slice()[slot], -1);
            assert_eq!(index.end.as_slice()[slot], -1);
            skipped += 1;
            continue;
        }
        assert_eq!(index.start.as_slice()[slot], cursor as i64);
        let charge = 2 + slot % 2;
        let mut expected = Vec::new();
        for ordinal in 1..view.len() {
            for z in 1..=charge {
                expected.push((87 + (view.len() - 1 - ordinal) * 3 + z - 1) as u8);
            }
            for z in 1..=charge {
                expected.push(((ordinal - 1) * 3 + z - 1) as u8);
            }
        }
        let end = cursor + expected.len();
        assert_eq!(index.end.as_slice()[slot], end as i64);
        assert_eq!(&arrays.annotation_id[cursor..end], expected);
        assert!(arrays.predicted_intensity.fp16().unwrap()[cursor..end]
            .iter()
            .all(|x| x.to_bits() == 0x3c00));
        cursor = end;
    }
    assert!(skipped > 0 && engine.calls > 3);
    assert_eq!(metrics.predicted_slots + skipped, options.row_count * 2);
    assert_eq!(metrics.fragments as usize, cursor);
    std::fs::remove_dir_all(path).unwrap();
}

#[test]
fn pipeline_failures_propagate_without_hanging() {
    let table = table();
    let tokenizer = Tokenizer::new(ModificationPolicy::Original, HashMap::new());
    let mut options = Options {
        row_start: 0,
        row_count: table.rows.len(),
        min_charge: 2,
        max_charge: 3,
        collision_energy: 30.0,
        batch_size: 1,
        top_k: 12,
        intensity_encoding: IntensityEncoding::Fp16,
        skip_unsupported: true,
    };
    let path = directory();
    let error = run(
        &table,
        &tokenizer,
        &mut CheckedEngine {
            calls: 0,
            fail: true,
        },
        &options,
        &path,
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("injected GPU failure"),
        "{error}"
    );
    std::fs::remove_dir_all(path).unwrap();
    options.skip_unsupported = false;
    let path = directory();
    let error = run(
        &table,
        &tokenizer,
        &mut CheckedEngine {
            calls: 0,
            fail: false,
        },
        &options,
        &path,
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("unsupported peptide row"),
        "{error}"
    );
    std::fs::remove_dir_all(path).unwrap();
}

#[test]
fn packed_writer_keeps_missing_empty_and_zero_predictions_distinct() {
    for encoding in [IntensityEncoding::Fp16, IntensityEncoding::Sqrt8Basepeak] {
        let path = directory();
        let mut writer = storage::Writer::new(&path, encoding).unwrap();
        let mut batch = Batch::new(3, 2, encoding);
        batch.slots = 4;
        batch.inputs = vec![Input { slot: 0 }, Input { slot: 2 }, Input { slot: 3 }];
        batch.output.resize(3);
        let stride = batch.output.stride();
        let width = batch.output.width();
        batch.output.data[0] = 2;
        batch.output.data[1 + 2 * width..3 + 2 * width].copy_from_slice(&[90, 0]);
        if width == 2 {
            batch.output.data[1..5].copy_from_slice(&[0, 0x48, 0, 0x80]); // 8, -0
        } else {
            batch.output.data[1..3].copy_from_slice(&[255, 0]);
        }
        // Row 1 is empty; row 2 has a present zero, distinct from missing slot 1.
        batch.output.data[2 * stride] = 1;
        batch.output.data[2 * stride + 1 + 2 * width] = 87;
        writer.append(&batch).unwrap();
        assert_eq!(writer.finish().unwrap(), (4, 3));
        let arrays = FragmentIntensityArrays::open(&path.join("arrays.mmappet")).unwrap();
        assert_eq!(arrays.annotation_id, [90, 0, 87]);
        if width == 2 {
            assert_eq!(
                arrays
                    .predicted_intensity
                    .fp16()
                    .unwrap()
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                [0x4800, 0x8000, 0]
            );
        } else {
            assert_eq!(arrays.predicted_intensity.codes().unwrap(), [255, 0, 0]);
        }
        let offsets = |column: &str| -> Vec<i64> {
            std::fs::read(path.join("index.mmappet").join(column))
                .unwrap()
                .chunks_exact(8)
                .map(|b| i64::from_le_bytes(b.try_into().unwrap()))
                .collect()
        };
        assert_eq!(offsets("0.bin"), [0, -1, 2, 2]);
        assert_eq!(offsets("1.bin"), [2, -1, 2, 3]);
        std::fs::remove_dir_all(path).unwrap();
    }
}

#[test]
fn packed_writer_rejects_invalid_counts_and_encoding() {
    for (count, output_encoding, reason) in [
        (255, IntensityEncoding::Fp16, "nonfinite"),
        (13, IntensityEncoding::Fp16, "fragment count"),
        (1, IntensityEncoding::Sqrt8Basepeak, "encoding mismatch"),
    ] {
        let path = directory();
        let mut writer = storage::Writer::new(&path, IntensityEncoding::Fp16).unwrap();
        let mut batch = Batch::new(1, 12, output_encoding);
        batch.slots = 1;
        batch.inputs.push(Input { slot: 0 });
        batch.output.resize(1);
        batch.output.data[0] = count;
        assert!(writer
            .append(&batch)
            .unwrap_err()
            .to_string()
            .contains(reason));
        drop(writer);
        std::fs::remove_dir_all(path).unwrap();
    }
    let mut predictions = engine::Predictions::new(usize::MAX, IntensityEncoding::Fp16);
    predictions.resize(1);
    assert_eq!(predictions.stride(), 1 + 174 * 3);
    assert!(predictions.row(1).is_err());
}
