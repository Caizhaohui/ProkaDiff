fn main() {
    let output = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output();
    let git_hash = match output {
        Ok(out) if out.status.success() => String::from_utf8_lossy(&out.stdout).trim().to_string(),
        _ => "UNKNOWN".to_string(),
    };
    println!("cargo:rustc-env=GIT_HASH={git_hash}");
}
