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
