use itertools::izip;
use itertools::Itertools;
use std::fs;
use std::process::Command;

#[test]
fn test_report() {
    assert!(Command::new("bash")
        .arg("-c")
        .arg("target/debug/fqc -q tests/resources/example.fastq > /tmp/report.html")
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());

    let result = fs::read_to_string("/tmp/report.html").unwrap();
    let expected = include_str!("expected/report.html");

    for (line, expected_line) in izip!(result.lines().sorted(), expected.lines().sorted()) {
        if !expected_line.contains("created")
            && !expected_line.contains("version")
            && !expected_line.contains("Spec")
        {
            assert_eq!(line, expected_line);
        }
    }
}

/// An empty FASTQ (a file with zero reads) must produce a valid, empty report
/// rather than panicking. This covers both a plain empty file (which needletail
/// reports as an `EmptyFile` error) and an empty gzipped file (which needletail
/// reports as an `Io` "failed to fill whole buffer" error); both previously
/// reached `parse_fastx_file(..).expect(..)` and aborted the process.
#[test]
fn test_empty_fastq_does_not_panic() {
    for input in [
        "tests/resources/empty.fastq",
        "tests/resources/empty.fastq.gz",
    ] {
        let summary =
            std::env::temp_dir().join(format!("fqc_empty_{}", input.rsplit('/').next().unwrap()));
        fs::create_dir_all(&summary).unwrap();

        let status = Command::new("target/debug/fqc")
            .args(["-q", input, "-s", summary.to_str().unwrap()])
            .stdout(std::process::Stdio::null())
            .status()
            .unwrap();

        assert!(status.success(), "fqc exited unsuccessfully on {}", input);
        assert!(
            summary.join("fastqc_data.txt").exists(),
            "fqc did not write fastqc_data.txt for {}",
            input
        );
    }
}

/// A FASTQ path that does not exist must surface an error instead of being
/// silently treated as an empty input. needletail reports a missing file as an
/// `Io` error, the same kind it uses for an empty gzipped file, so this guards
/// the `is_file` distinction `process` relies on to tell the two apart.
#[test]
fn test_missing_fastq_errors() {
    let missing = std::env::temp_dir().join("fqc_nonexistent_input.fastq.gz");
    let _ = fs::remove_file(&missing);

    let status = Command::new("target/debug/fqc")
        .args(["-q", missing.to_str().unwrap()])
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();

    assert!(
        !status.success(),
        "fqc should fail on a missing input file, but it succeeded"
    );
}

/// An empty gzipped *stream* must be handled gracefully too. Process
/// substitution exposes the input as a non-regular file under `/dev/fd`, so a
/// `std::path::Path::is_file` check (false for pipes, FIFOs, and character
/// devices such as `/dev/stdin`) would wrongly reject a valid empty stream.
#[test]
fn test_empty_gzip_stream_does_not_panic() {
    let status = Command::new("bash")
        .arg("-c")
        .arg("target/debug/fqc -q <(printf '' | gzip) > /dev/null")
        .status()
        .unwrap();

    assert!(status.success(), "fqc failed on an empty gzipped stream");
}
