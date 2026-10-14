//! This process's CPU time (Linux only), for how busy a phase kept the cores.

/// CPU seconds this process has used, user and system, over all its threads; `None`
/// where `/proc` is missing.
pub fn cpu_seconds() -> Option<f64> {
    // `/proc/<pid>/stat` counts in USER_HZ ticks, 100 per second on Linux.
    const TICKS_PER_SECOND: f64 = 100.0;
    let stat = std::fs::read_to_string("/proc/self/stat").ok()?;
    // Fields after the parenthesised command name, which may itself contain spaces.
    let fields: Vec<&str> = stat.rsplit_once(')')?.1.split_whitespace().collect();
    let user: f64 = fields.get(11)?.parse().ok()?;
    let system: f64 = fields.get(12)?.parse().ok()?;
    Some((user + system) / TICKS_PER_SECOND)
}

/// Share of all cores a phase kept busy, from its CPU seconds and wall seconds.
pub fn utilisation(cpu_seconds: f64, wall_seconds: f64) -> f64 {
    cpu_seconds / (wall_seconds * rayon::current_num_threads() as f64).max(f64::EPSILON)
}
