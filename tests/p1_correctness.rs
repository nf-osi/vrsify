//! Regression coverage for the pre-publication review findings (P1 and P2).

use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::Value;

const SEQUENCE: &str = "CGTAAAAAACGTACGT";
const VCF_HEADER: &str = "##fileformat=VCFv4.2\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
const MAF_HEADER: &str = "Chromosome\tStart_Position\tReference_Allele\t\
    Tumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n";

struct Fixture(PathBuf);

impl Fixture {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let dir = std::env::temp_dir().join(format!(
            "vrsify_p1_{}_{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        std::fs::write(dir.join("ref.fa"), format!(">chr1\n{SEQUENCE}\n")).unwrap();
        std::fs::write(
            dir.join("map.tsv"),
            format!("chr1\t{}\n", vrsify::refget::refget_accession(SEQUENCE.as_bytes())),
        )
        .unwrap();
        Self(dir)
    }

    fn run(&self, maf: bool, rows: &str, reference: Option<&str>, extra: &[&str]) -> Output {
        let input = self.0.join(if maf { "input.maf" } else { "input.vcf" });
        let header = if maf { MAF_HEADER } else { VCF_HEADER };
        std::fs::write(&input, format!("{header}{rows}")).unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_vrsify"));
        if maf {
            cmd.arg("maf");
        }
        cmd.arg(if maf { "--maf" } else { "--vcf" })
            .arg(input)
            .arg("--seqmap")
            .arg(self.0.join("map.tsv"))
            .arg("--out-alleles")
            .arg(self.0.join("alleles.ndjson"))
            .arg("--out-observations")
            .arg(self.0.join("obs.ndjson"))
            .args(extra);
        if let Some(reference) = reference {
            cmd.arg("--reference").arg(self.0.join(reference));
        }
        cmd.output().unwrap()
    }

    fn json(&self, name: &str) -> Vec<Value> {
        std::fs::read_to_string(self.0.join(name))
            .unwrap()
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect()
    }

    fn assert_rejected(&self, output: Output, message: &str) {
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "expected an error, not a panic: {stderr}");
        assert!(stderr.contains(message), "{stderr}");
        assert!(self.json("alleles.ndjson").is_empty(), "invalid allele was emitted");
        assert!(self.json("obs.ndjson").is_empty(), "invalid observation was emitted");
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

fn row(maf: bool, pos: u64, reference: &str, alt: &str) -> String {
    if maf {
        format!("chr1\t{pos}\t{reference}\t{reference}\t{alt}\tS1\n")
    } else {
        format!("chr1\t{pos}\t.\t{reference}\t{alt}\t.\tPASS\t.\tGT\t0/1\n")
    }
}

#[test]
fn missing_maf_alleles_never_become_empty_alleles() {
    let f = Fixture::new();
    for missing in ["", "NA", ".", " "] {
        for extra in [&[][..], &["--strict"][..]] {
            for reference in [None, Some("ref.fa")] {
                let output = f.run(true, &row(true, 4, "A", missing), reference, extra);
                f.assert_rejected(output, "MAF row 1: missing Tumor_Seq_Allele2");
                let output = f.run(true, &row(true, 4, missing, "T"), reference, extra);
                f.assert_rejected(output, "MAF row 1: missing Reference_Allele");
            }
        }
    }
}

#[test]
fn explicit_empty_maf_alleles_remain_valid() {
    let f = Fixture::new();
    for (reference, alt) in [("A", "-"), ("-", "A")] {
        let output = f.run(true, &row(true, 4, reference, alt), Some("ref.fa"), &[]);
        assert!(output.status.success(), "{:?}", output);
        let alleles = f.json("alleles.ndjson");
        assert_eq!(alleles.len(), 1);
        assert!(alleles[0].get("fullyJustified").is_none());
        assert_eq!(f.json("obs.ndjson")[0]["alternateBases"], alt);
    }
}

#[test]
fn ref_mismatches_fail_in_both_formats_including_former_panic() {
    let f = Fixture::new();
    for maf in [false, true] {
        // FASTA has A here. C>A formerly entered normalization with an empty unit.
        for alt in ["G", "A"] {
            let output = f.run(maf, &row(maf, 4, "C", alt), Some("ref.fa"), &[]);
            f.assert_rejected(output, "REF mismatch");
        }
    }
}

#[test]
fn fasta_must_match_seqmap_even_when_local_ref_bases_match() {
    let f = Fixture::new();
    // Only the LAST base differs; the input variant's REF still agrees locally.
    std::fs::write(f.0.join("wrong.fa"), ">chr1\nCGTAAAAAACGTACGA\n").unwrap();
    for maf in [false, true] {
        let output = f.run(maf, &row(maf, 4, "A", "T"), Some("wrong.fa"), &[]);
        f.assert_rejected(output, "does not match seqmap accession");
    }
}

#[test]
fn reference_bounds_are_checked_before_normalization() {
    let f = Fixture::new();
    for maf in [false, true] {
        for (pos, reference, alt) in [(17, "A", "T"), (16, "TA", "T"), (i64::MAX as u64, "A", "T")] {
            let output = f.run(maf, &row(maf, pos, reference, alt), Some("ref.fa"), &[]);
            f.assert_rejected(output, "outside reference bounds");
        }
    }
    let output = f.run(true, &row(true, 17, "-", "A"), Some("ref.fa"), &[]);
    f.assert_rejected(output, "outside reference bounds");
}

#[test]
fn supplied_reference_cannot_silently_fall_back_on_missing_contig() {
    let f = Fixture::new();
    std::fs::write(f.0.join("other.fa"), ">other\nACGT\n").unwrap();
    for maf in [false, true] {
        let output = f.run(maf, &row(maf, 4, "A", "T"), Some("other.fa"), &[]);
        f.assert_rejected(output, "reference FASTA is missing contig 'chr1'");
    }
}

#[test]
fn vcf_indels_without_reference_are_flagged_and_warned() {
    let f = Fixture::new();
    let rows = row(false, 4, "AA", "A") + &row(false, 4, "A", "AA") + &row(false, 4, "A", "T");
    let output = f.run(false, &rows, None, &[]);
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("2 indel allele(s) are NOT fully justified"));
    let alleles = f.json("alleles.ndjson");
    assert_eq!(alleles.len(), 3);
    assert_eq!(alleles[0]["fullyJustified"], false);
    assert_eq!(alleles[1]["fullyJustified"], false);
    assert!(alleles[2].get("fullyJustified").is_none());
    let observations = f.json("obs.ndjson");
    for (allele, obs) in alleles.iter().zip(&observations) {
        assert_eq!(allele["id"], obs["variant"]);
        let model: vrsify::vrs::Allele = serde_json::from_value(allele.clone()).unwrap();
        assert_eq!(allele["id"], model.ga4gh_id(), "marker must not affect digest");
    }

    let output = f.run(false, &rows, Some("ref.fa"), &[]);
    assert!(output.status.success());
    assert!(!String::from_utf8_lossy(&output.stderr).contains("NOT fully justified"));
    let normalized = f.json("alleles.ndjson");
    assert!(normalized.iter().all(|a| a.get("fullyJustified").is_none()));
    assert_eq!(normalized[0]["id"], "ga4gh:VA.F0wcxqL5qakH_o89aqX19pUYoBRVDqmI");
    assert_ne!(alleles[0]["id"], normalized[0]["id"]);
}

#[test]
fn reference_validation_accepts_case_insensitive_ref_bases() {
    let f = Fixture::new();
    for maf in [false, true] {
        let output = f.run(maf, &row(maf, 4, "a", "t"), Some("ref.fa"), &[]);
        assert!(output.status.success(), "{:?}", output);
        assert_eq!(f.json("alleles.ndjson")[0]["state"]["sequence"], "T");
    }
}

/// P2: the seqmap is keyed `chr1`, the source file writes `1`. The MAF path resolved
/// that alias; the VCF path did an exact lookup and quietly emitted nothing.
#[test]
fn contig_aliases_resolve_in_both_formats() {
    let f = Fixture::new();
    for maf in [false, true] {
        for reference in [None, Some("ref.fa")] {
            let row = if maf {
                "1\t4\tA\tA\tT\tS1\n".to_string()
            } else {
                "1\t4\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n".to_string()
            };
            let output = f.run(maf, &row, reference, &[]);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(output.status.success(), "{stderr}");
            assert!(!stderr.contains("missing from seqmap"), "{stderr}");
            // Same id as the `chr1`-spelled row: the alias is a naming difference, not
            // a different variant.
            let alleles = f.json("alleles.ndjson");
            assert_eq!(alleles.len(), 1, "alias contig produced no allele: {stderr}");
            assert_eq!(f.json("obs.ndjson")[0]["variant"], alleles[0]["id"]);
            assert_eq!(alleles[0]["location"]["start"], 3);
        }
    }
}

/// P2: VCF nucleotide strings are case-insensitive, so the state must be uppercased
/// before hashing — otherwise `A>t` and `A>T` mint two ids for one variant. With a
/// `--reference` this came free from `normalize`; without one it did not.
#[test]
fn lowercase_vcf_alleles_get_the_uppercase_id() {
    let f = Fixture::new();
    let rows = row(false, 4, "A", "T") + &row(false, 4, "a", "t");
    let output = f.run(false, &rows, None, &[]);
    assert!(output.status.success(), "{:?}", output);
    let alleles = f.json("alleles.ndjson");
    assert_eq!(alleles.len(), 1, "case variants must collapse: {alleles:#?}");
    assert_eq!(alleles[0]["state"]["sequence"], "T");
    let observations = f.json("obs.ndjson");
    assert_eq!(observations.len(), 2);
    for obs in &observations {
        assert_eq!(obs["variant"], alleles[0]["id"]);
    }
    // The observation keeps the row as written; only the VRS state is normalized.
    assert_eq!(observations[1]["alternateBases"], "t");
}

/// P2: `GT=1/.` is a diploid call with one uncalled allele — it could be `1/0` or
/// `1/1`. Counting only the called alleles made it look haploid, i.e. `hemizygous`.
#[test]
fn partially_missing_genotype_is_not_hemizygous() {
    let f = Fixture::new();
    let gt_row = |gt: &str| format!("chr1\t4\t.\tA\tT\t.\tPASS\t.\tGT\t{gt}\n");
    for gt in ["1/.", "./1", ".|1", "1|."] {
        let output = f.run(false, &gt_row(gt), None, &[]);
        assert!(output.status.success(), "{:?}", output);
        assert_eq!(f.json("obs.ndjson")[0]["zygosity"], "unknown", "GT={gt}");
    }
    // Genuinely haploid and fully-called genotypes are unaffected.
    for (gt, zygosity) in [("1", "hemizygous"), ("1/1", "homozygous"), ("0/1", "heterozygous")] {
        let output = f.run(false, &gt_row(gt), None, &[]);
        assert!(output.status.success(), "{:?}", output);
        assert_eq!(f.json("obs.ndjson")[0]["zygosity"], zygosity, "GT={gt}");
    }
}

/// P2 follow-up: a rejected row used `unknown` in its local variant id when the MAF
/// omitted NCBI_Build, even though the resolved seqmap supplied an assembly that was
/// copied onto the observation. The node and observation must derive assembly once.
#[test]
fn rejected_variant_and_observation_share_seqmap_assembly() {
    let f = Fixture::new();
    std::fs::write(
        f.0.join("map.tsv"),
        format!(
            "chr1\t{}\tSYN1\t1\n",
            vrsify::refget::refget_accession(SEQUENCE.as_bytes())
        ),
    )
    .unwrap();

    let output = f.run(true, &row(true, 4, "A", "Z"), None, &[]);
    assert!(output.status.success(), "{:?}", output);

    let alleles = f.json("alleles.ndjson");
    let observations = f.json("obs.ndjson");
    assert_eq!(alleles.len(), 1);
    assert_eq!(observations.len(), 1);
    assert_eq!(alleles[0]["id"], "nf:variant/SYN1:chr1:4:A:Z");
    assert_eq!(alleles[0]["assemblyId"], "SYN1");
    assert_eq!(observations[0]["variant"], alleles[0]["id"]);
    assert_eq!(observations[0]["assemblyId"], alleles[0]["assemblyId"]);
}
