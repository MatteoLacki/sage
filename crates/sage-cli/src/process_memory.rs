//! This process's peak resident memory and the machine's available memory (Linux only).

/// Writing "5" to `/proc/self/clear_refs` resets `VmHWM` to the current RSS.
pub fn reset_peak_rss() -> std::io::Result<()> {
    std::fs::write("/proc/self/clear_refs", "5")
}

pub fn peak_rss_bytes() -> anyhow::Result<usize> {
    kilobytes_field("/proc/self/status", "VmHWM:").map(|kilobytes| kilobytes * 1024)
}

/// Memory the kernel estimates can be allocated without swapping (`MemAvailable`).
pub fn available_memory_bytes() -> anyhow::Result<usize> {
    kilobytes_field("/proc/meminfo", "MemAvailable:").map(|kilobytes| kilobytes * 1024)
}

fn kilobytes_field(path: &str, field: &str) -> anyhow::Result<usize> {
    let contents = std::fs::read_to_string(path)?;
    let line = contents
        .lines()
        .find(|line| line.starts_with(field))
        .ok_or_else(|| anyhow::anyhow!("no {field} line in {path}"))?;
    let kilobytes: usize = line
        .trim_start_matches(field)
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()?;
    Ok(kilobytes)
}
