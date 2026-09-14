from pathlib import Path

root = Path("target/simd-benchmark").resolve()
p = Path("crates/sage/src/database.rs")
s = p.read_text()
assert "report_page_scan_statistics" not in s, "Already instrumented"
(root / "database.uninstrumented.rs").write_text(s)
s = s.replace(
    "                page_window,\n            }",
    "                page_window,\n                histogram,\n                lanes,\n            }",
)
s = s.replace(
    "                i = j;\n",
    """                let windows_on_page = (j - i) as u64;
                histogram[masses.len().min(64)] += windows_on_page;
                lanes[0] += masses.len() as u64 * windows_on_page;
                lanes[1] += (masses.len() / 8 * 8) as u64 * windows_on_page;
                i = j;
""",
)
s = s.replace(
    "    page_window: Vec<(usize, usize)>,\n}",
    "    page_window: Vec<(usize, usize)>,\n    histogram: [u64; 65],\n    lanes: [u64; 2],\n}",
)
s = s.replace(
    "#[derive(Default)]\nstruct PageSearchScratch", "struct PageSearchScratch"
)
s = s.replace(
    "thread_local! {\n    static PAGE_SEARCH_SCRATCH",
    """impl Default for PageSearchScratch {
    fn default() -> Self { Self { bounds: Vec::new(), page_window: Vec::new(), histogram: [0; 65], lanes: [0; 2] } }
}

pub fn report_page_scan_statistics() {
    let workers = rayon::broadcast(|_| PAGE_SEARCH_SCRATCH.with(|s| {
        let s = s.borrow();
        (s.histogram, s.lanes)
    }));
    let mut histogram = [0u64; 65];
    let mut lanes = [0u64; 2];
    for (counts, scanned) in workers {
        for (out, value) in histogram.iter_mut().zip(counts) { *out += value; }
        for (out, value) in lanes.iter_mut().zip(scanned) { *out += value; }
    }
    log::info!("PAGE_SCAN_HISTOGRAM {:?}; LANES {:?}", histogram, lanes);
}

thread_local! {
    static PAGE_SEARCH_SCRATCH""",
)
p.write_text(s)
p = Path("crates/sage-cli/src/runner.rs")
s = p.read_text()
(root / "runner.uninstrumented.rs").write_text(s)
s = s.replace(
    "        features\n    }\n\n    fn complete_features",
    "        sage_core::database::report_page_scan_statistics();\n        features\n    }\n\n    fn complete_features",
)
p.write_text(s)
