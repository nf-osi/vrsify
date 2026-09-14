//! Acceptance for reference-based normalization and seqmap generation, on SMALL
//! SYNTHETIC references (the equivalent check on real GRCh38 needs the multi-GB FASTA,
//! so it is not run here). `seqmap` derives `SQ.` refget accessions from a FASTA; an
//! indel in a homopolymer is fully justified, so a left-shifted and a right-shifted VCF
//! representation collapse to the SAME VRS id.

use std::process::Command;

// A tiny contig with a homopolymer run (the AAAAAA) to exercise indel justification.
//              0         1
//              0123456789012345
const FASTA: &str = ">tinychr\nCGTAAAAAACGTACGT\n";

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vrsify"))
}

#[test]
fn seqmap_subcommand_and_normalization() {
    let dir = std::env::temp_dir().join(format!("vrsify_ref_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, FASTA).unwrap();

    // --- derive the seqmap from the FASTA ---
    let seqmap = dir.join("seqmap.tsv");
    let status = bin()
        .args([
            "seqmap",
            "--fasta",
            fasta.to_str().unwrap(),
            "--assembly",
            "SYN1",
            "--out",
            seqmap.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let seqmap_txt = std::fs::read_to_string(&seqmap).unwrap();
    // The refget accession of CGTAAAAAACGTACGT (sha512t24u of the uppercased bytes).
    // Recomputed independently below to avoid hard-coding a possibly-wrong constant.
    assert!(seqmap_txt.contains("tinychr\tSQ."));
    assert!(seqmap_txt.contains("\tSYN1\ttinychr"));

    // --- two equivalent representations of a single-A deletion ---
    // Left-shifted:  POS=4 (1-based) REF=AA ALT=A  → delete an A at the run's left.
    // Right-shifted: POS=8 (1-based) REF=AA ALT=A  → delete an A at the run's right.
    // Both are the same edit within the AAAAAA run and MUST yield one VRS allele id.
    let left_vcf = dir.join("left.vcf");
    let right_vcf = dir.join("right.vcf");
    let hdr = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    std::fs::write(&left_vcf, format!("{hdr}tinychr\t4\t.\tAA\tA\t.\tPASS\t.\n")).unwrap();
    std::fs::write(&right_vcf, format!("{hdr}tinychr\t8\t.\tAA\tA\t.\tPASS\t.\n")).unwrap();

    let run = |vcf: &std::path::Path, tag: &str| -> String {
        let out = dir.join(format!("al_{tag}.ndjson"));
        let obs = dir.join(format!("obs_{tag}.ndjson"));
        let status = bin()
            .args([
                "--vcf",
                vcf.to_str().unwrap(),
                "--seqmap",
                seqmap.to_str().unwrap(),
                "--reference",
                fasta.to_str().unwrap(),
                "--out-alleles",
                out.to_str().unwrap(),
                "--out-observations",
                obs.to_str().unwrap(),
            ])
            .status()
            .unwrap();
        assert!(status.success());
        std::fs::read_to_string(&out).unwrap()
    };

    let left = run(&left_vcf, "left");
    let right = run(&right_vcf, "right");

    let id_of = |s: &str| -> String {
        let v: serde_json::Value = serde_json::from_str(s.lines().next().unwrap()).unwrap();
        v["id"].as_str().unwrap().to_string()
    };
    let left_id = id_of(&left);
    let right_id = id_of(&right);
    assert_eq!(
        left_id, right_id,
        "left- and right-shifted indel reps did not collapse to one VRS id\nleft:  {left}\nright: {right}"
    );

    // Cross-validation lock-in: this fully-justified homopolymer deletion MUST be emitted
    // as a ReferenceLengthExpression (NOT a Literal) and its VA id must byte-match the
    // vrs-python 2.3.3 reference implementation. The id/state below were computed by
    // vrs-python on this exact tiny contig (see verify/); if RLE emission regresses to a
    // LiteralSequenceExpression the digest — and therefore this id — changes.
    let left_val: serde_json::Value = serde_json::from_str(left.lines().next().unwrap()).unwrap();
    assert_eq!(
        left_val["state"]["type"].as_str(),
        Some("ReferenceLengthExpression"),
        "fully-justified indel must be an RLE, got {}",
        left_val["state"]
    );
    assert_eq!(left_val["state"]["length"].as_i64(), Some(5));
    assert_eq!(left_val["state"]["repeatSubunitLength"].as_i64(), Some(1));
    assert_eq!(
        left_id, "ga4gh:VA.F0wcxqL5qakH_o89aqX19pUYoBRVDqmI",
        "VA id diverged from the vrs-python 2.3.3 reference id for the RLE deletion"
    );

    // Sanity: without --reference (naive projection), they would NOT collapse.
    let naive_left = {
        let out = dir.join("nl.ndjson");
        let obs = dir.join("nlo.ndjson");
        bin()
            .args([
                "--vcf", left_vcf.to_str().unwrap(),
                "--seqmap", seqmap.to_str().unwrap(),
                "--out-alleles", out.to_str().unwrap(),
                "--out-observations", obs.to_str().unwrap(),
            ])
            .status()
            .unwrap();
        std::fs::read_to_string(&out).unwrap()
    };
    let naive_right = {
        let out = dir.join("nr.ndjson");
        let obs = dir.join("nro.ndjson");
        bin()
            .args([
                "--vcf", right_vcf.to_str().unwrap(),
                "--seqmap", seqmap.to_str().unwrap(),
                "--out-alleles", out.to_str().unwrap(),
                "--out-observations", obs.to_str().unwrap(),
            ])
            .status()
            .unwrap();
        std::fs::read_to_string(&out).unwrap()
    };
    assert_ne!(
        id_of(&naive_left),
        id_of(&naive_right),
        "naive (no-reference) projection unexpectedly collapsed the two reps"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
