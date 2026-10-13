//! This process's peak resident memory, for the sizing and verification binaries (Linux only).

/// Writing "5" to `/proc/self/clear_refs` resets `VmHWM` to the current RSS.
pub fn reset_peak_rss() -> std::io::Result<()> {
    std::fs::write("/proc/self/clear_refs", "5")
}

pub fn peak_rss_bytes() -> anyhow::Result<usize> {
    let status = std::fs::read_to_string("/proc/self/status")?;
    let line = status
        .lines()
        .find(|line| line.starts_with("VmHWM:"))
        .ok_or_else(|| anyhow::anyhow!("no VmHWM line in /proc/self/status"))?;
    let kilobytes: usize = line
        .trim_start_matches("VmHWM:")
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()?;
    Ok(kilobytes * 1024)
}
