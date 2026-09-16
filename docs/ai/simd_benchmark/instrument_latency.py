"""Instrument `page_search_batch` to separate per-page peptide-ID bound
resolution from the fragment-mass scan that follows it.

Counts page groups and windows and accumulates TSC cycles in each phase, so
the two can be weighed against the search-phase wall time. Timing brackets sit
around whole phases rather than individual pages, to keep the serializing
`lfence` pairs off the hot per-page path.

Saves the original sources; `build_latency_diagnostic.py` restores them.
"""

from pathlib import Path

root = Path("target/simd-benchmark").resolve()

p = Path("crates/sage/src/database.rs")
s = p.read_text()
assert "report_page_latency_statistics" not in s, "Already instrumented"
(root / "database.unlatency.rs").write_text(s)

s = s.replace(
    "                page_window,\n                groups,\n                page_bounds,\n            }",
    "                page_window,\n                groups,\n                page_bounds,\n                stats,\n            }",
)
s = s.replace(
    "    /// `(inner_left, inner_right)` per entry in `groups`.\n    page_bounds: Vec<(u32, u32)>,\n}",
    "    /// `(inner_left, inner_right)` per entry in `groups`.\n    page_bounds: Vec<(u32, u32)>,\n    stats: [u64; 4],\n}",
)

old_resolve = """            resolve_page_peptide_bounds(
                &self.db.fragment_peptide_indices,
                self.db.bucket_size,
                self.db.size(),
                groups,
                exact_pre_idx_lo,
                self.pre_idx_hi,
                page_bounds,
            );
"""
new_resolve = """            let t0 = unsafe { core::arch::x86_64::_mm_lfence(); core::arch::x86_64::_rdtsc() };
            resolve_page_peptide_bounds(
                &self.db.fragment_peptide_indices,
                self.db.bucket_size,
                self.db.size(),
                groups,
                exact_pre_idx_lo,
                self.pre_idx_hi,
                page_bounds,
            );
            let t1 = unsafe { core::arch::x86_64::_mm_lfence(); core::arch::x86_64::_rdtsc() };
            stats[0] += groups.len() as u64;
            stats[1] += page_window.len() as u64;
            stats[2] += t1 - t0;
"""
assert old_resolve in s, "resolve call not found"
s = s.replace(old_resolve, new_resolve)

old_tail = """                    scan_masses_scalar(masses, fragment_lo, fragment_hi, &mut emit);
                    j += 1;
                }
            }
        });
"""
new_tail = """                    scan_masses_scalar(masses, fragment_lo, fragment_hi, &mut emit);
                    j += 1;
                }
            }
            let t2 = unsafe { core::arch::x86_64::_mm_lfence(); core::arch::x86_64::_rdtsc() };
            stats[3] += t2 - t1;
        });
"""
assert old_tail in s, "scan tail not found"
s = s.replace(old_tail, new_tail)

s = s.replace("#[derive(Default)]\nstruct PageSearchScratch", "struct PageSearchScratch")
s = s.replace(
    "thread_local! {\n    static PAGE_SEARCH_SCRATCH",
    """impl Default for PageSearchScratch {
    fn default() -> Self {
        Self {
            bounds: Vec::new(),
            page_window: Vec::new(),
            groups: Vec::new(),
            page_bounds: Vec::new(),
            stats: [0; 4],
        }
    }
}

pub fn report_page_latency_statistics() {
    let workers = rayon::broadcast(|_| PAGE_SEARCH_SCRATCH.with(|s| s.borrow().stats));
    let mut stats = [0u64; 4];
    for worker in workers {
        for (out, value) in stats.iter_mut().zip(worker) {
            *out += value;
        }
    }
    log::info!(
        "PAGE_LATENCY page_groups={} windows={} bound_cycles={} scan_cycles={}",
        stats[0], stats[1], stats[2], stats[3]
    );
}

thread_local! {
    static PAGE_SEARCH_SCRATCH""",
)
p.write_text(s)

p = Path("crates/sage-cli/src/runner.rs")
s = p.read_text()
(root / "runner.unlatency.rs").write_text(s)
old = "        features\n    }\n\n    fn complete_features"
assert old in s, "runner hook site not found"
s = s.replace(
    old,
    "        sage_core::database::report_page_latency_statistics();\n"
    "        features\n    }\n\n    fn complete_features",
)
p.write_text(s)
print("instrumented")
