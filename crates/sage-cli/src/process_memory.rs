//! This process's resident memory and the machine's available memory (Linux only).

/// Writing "5" to `/proc/self/clear_refs` resets `VmHWM` to the current RSS.
pub fn reset_peak_rss() -> std::io::Result<()> {
    std::fs::write("/proc/self/clear_refs", "5")
}

pub fn peak_rss_bytes() -> anyhow::Result<usize> {
    kilobytes_field("/proc/self/status", "VmHWM:").map(|kilobytes| kilobytes * 1024)
}

/// Private resident memory (`RssAnon`): the heap and other pages the kernel can't drop,
/// unlike mapped files, which `peak_rss_bytes` also counts.
pub fn anonymous_rss_bytes() -> anyhow::Result<usize> {
    kilobytes_field("/proc/self/status", "RssAnon:").map(|kilobytes| kilobytes * 1024)
}

/// Logs this process's anonymous and peak resident memory; silent where `/proc` is missing.
pub fn log_memory(checkpoint: &str) {
    if let (Ok(anonymous), Ok(peak)) = (anonymous_rss_bytes(), peak_rss_bytes()) {
        log::info!(
            "- memory, {checkpoint}: anonymous {:.2} GB, peak RSS {:.2} GB",
            anonymous as f64 / 1e9,
            peak as f64 / 1e9
        );
    }
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
