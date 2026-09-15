//! Acceptance: MAF (cBioPortal `data_mutations.txt`) → VRS alleles + somatic
//! observations.
//!
//! The contract this locks in is **cross-format identity**: the same variant, whether it
//! arrives as a VCF record or a MAF row, must get the *same* `ga4gh:VA.` id. That is the
//! whole reason for minting VRS ids in the KG (issue nf-osi/kg-pipeline#95) — otherwise
//! a variant ingested from cBioPortal would not join to one ingested from a VCF, to
//! ClinVar, or to gnomAD.
//!
//! * `maf_snv_reaches_the_ga4gh_golden_id` — a MAF SNP row for rs7412 produces the
//!   GA4GH `vrs@2.0` golden fixture id, byte-for-byte, the same one the VCF path
//!   produces from `examples/sample.vcf`.
//! * `maf_and_vcf_indel_agree_after_normalization` — a deletion and an insertion,
//!   expressed in MAF's trimmed `-` form and in VCF's anchored form, collapse to one id
//!   (the vrs-python 2.3.3 cross-validated id) once `--reference` is supplied.
//! * `indel_without_reference_is_flagged` — and when it is *not* supplied, the id
//!   diverges and the allele says so rather than lying.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The same tiny homopolymer contig `tests/reference_m3_m4.rs` uses, so the indel ids
/// below are comparable with the VCF-path cross-validation there.
//                        0         1
//                        0123456789012345
const FASTA: &str = ">tinychr\nCGTAAAAAACGTACGT\n";

/// vrs-python 2.3.3's VA id for the fully-justified single-A deletion in that run.
const GOLDEN_DEL_ID: &str = "ga4gh:VA.F0wcxqL5qakH_o89aqX19pUYoBRVDqmI";

/// GA4GH `vrs@2.0` `models.yaml` golden id for rs7412 (GRCh38 chr19:44908822 C>T).
const GOLDEN_RS7412_ID: &str = "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt";

const MAF_HEADER: &str = "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\t\
    Start_Position\tEnd_Position\tStrand\tConsequence\tVariant_Classification\tVariant_Type\t\
    Reference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tTumor_Sample_Barcode\t\
    Matched_Norm_Sample_Barcode\tHGVSc\tHGVSp_Short\tt_depth\tt_ref_count\tt_alt_count";

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vrsify"))
}

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("vrsify_maf_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Run `vrsify maf` and return (alleles NDJSON, observations NDJSON, stderr).
fn run_maf(dir: &Path, maf: &Path, seqmap: &Path, extra: &[&str]) -> (String, String, String) {
    let alleles = dir.join("alleles.ndjson");
    let obs = dir.join("obs.ndjson");
    let out = bin()
        .args(["maf", "--maf"])
        .arg(maf)
        .arg("--seqmap")
        .arg(seqmap)
        .arg("--out-alleles")
        .arg(&alleles)
        .arg("--out-observations")
        .arg(&obs)
        .args(extra)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(out.status.success(), "vrsify maf failed: {stderr}");
    (
        std::fs::read_to_string(&alleles).unwrap(),
        std::fs::read_to_string(&obs).unwrap(),
        stderr,
    )
}

fn json_lines(s: &str) -> Vec<serde_json::Value> {
    s.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("NDJSON line parses"))
        .collect()
}

fn id_of(value: &serde_json::Value) -> &str {
    value["id"].as_str().expect("allele has an id")
}

/// Write a seqmap for the synthetic contig by running the `seqmap` subcommand, so the
/// refget accession is derived the same way a real run derives it.
fn synthetic_seqmap(dir: &Path) -> (PathBuf, PathBuf) {
    let fasta = dir.join("tiny.fa");
    std::fs::write(&fasta, FASTA).unwrap();
    let seqmap = dir.join("seqmap.tsv");
    let status = bin()
        .args(["seqmap", "--fasta"])
        .arg(&fasta)
        .arg("--assembly")
        .arg("SYN1")
        .arg("--out")
        .arg(&seqmap)
        .status()
        .unwrap();
    assert!(status.success());
    (fasta, seqmap)
}

#[test]
fn maf_snv_reaches_the_ga4gh_golden_id() {
    let dir = scratch("snv");
    // Two tumor samples carry the same rs7412 allele — the recurrence the variant layer
    // exists to answer ("has this variant been seen in a case?").
    //
    // Note `Chromosome` is `19` while `examples/seqmap.tsv` is keyed `chr19`: cBioPortal
    // MAFs drop the `chr` prefix a UCSC-derived seqmap keeps, and the lookup bridges it.
    let maf = dir.join("data_mutations.txt");
    std::fs::write(
        &maf,
        format!(
            "#genome_nexus_version: 1.0.2\n#isoform: mskcc\n{MAF_HEADER}\n\
             APOE\t348\tSage\tGRCh38\t19\t44908822\t44908822\t+\tmissense_variant\t\
             Missense_Mutation\tSNP\tC\tC\tT\trs7412\tJH-2-001-A\tJH-2-001-N\t\
             ENST00000252486.9:c.388T>C\tp.C130R\t60\t40\t20\n\
             APOE\t348\tSage\tGRCh38\t19\t44908822\t44908822\t+\tmissense_variant\t\
             Missense_Mutation\tSNP\tC\tC\tT\trs7412\tJH-2-002-A\tJH-2-002-N\t\
             ENST00000252486.9:c.388T>C\tp.C130R\t\t\t\n"
        ),
    )
    .unwrap();

    let seqmap = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/seqmap.tsv");
    // No --reference needed: the trimmed projection is already exact for substitutions.
    let (alleles, obs, _) = run_maf(
        &dir,
        &maf,
        &seqmap,
        &["--study-id", "nst_nfosi_ntap", "--source", "syn12345/maf"],
    );

    let alleles = json_lines(&alleles);
    assert_eq!(
        alleles.len(),
        1,
        "the recurrent allele must be emitted once, got {alleles:#?}"
    );
    assert_eq!(
        id_of(&alleles[0]),
        GOLDEN_RS7412_ID,
        "MAF path diverged from the GA4GH golden id the VCF path produces"
    );
    // Nothing sample-specific may leak onto the context-free allele.
    let allele_keys: Vec<&String> = alleles[0].as_object().unwrap().keys().collect();
    assert_eq!(allele_keys, ["id", "location", "state", "type"]);

    let obs = json_lines(&obs);
    assert_eq!(obs.len(), 2, "one observation per tumor sample");
    for o in &obs {
        assert_eq!(o["type"], "VariantObservation");
        assert_eq!(o["variant"], GOLDEN_RS7412_ID);
        assert_eq!(o["assemblyId"], "GRCh38");
        assert_eq!(o["referenceName"], "19", "denormalized Beacon referenceName");
        assert_eq!(o["studyId"], "nst_nfosi_ntap");
        // NF-OSI MAFs leave Mutation_Status blank; the profile is somatic by construction.
        assert_eq!(o["mutationStatus"], "Somatic");
        assert_eq!(o["affectedGeneSymbol"], "APOE");
        assert_eq!(o["entrezGeneId"], "348");
        assert_eq!(o["aminoacidChange"], "p.C130R", "HGVSp_Short is searchable");
        assert_eq!(o["variantClassification"], "Missense_Mutation");
        // Consequence is a list so each SO term can be mapped through SSSOM.
        assert_eq!(o["molecularConsequence"], serde_json::json!(["missense_variant"]));
        assert_eq!(o["dbsnpId"], serde_json::json!(["rs7412"]));
    }
    assert_eq!(obs[0]["tumorSampleBarcode"], "JH-2-001-A");
    assert_eq!(obs[0]["matchedNormalSampleBarcode"], "JH-2-001-N");
    assert_eq!(obs[0]["tumorAltCount"], 20);
    assert_eq!(obs[0]["variantAlleleFrequency"], 20.0 / 60.0);
    // Second row has no depth columns filled — absent, not zero.
    assert_eq!(obs[1]["tumorSampleBarcode"], "JH-2-002-A");
    assert!(obs[1].get("variantAlleleFrequency").is_none());
    assert!(obs[1].get("tumorAltCount").is_none());

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn maf_and_vcf_indel_agree_after_normalization() {
    let dir = scratch("indel");
    let (fasta, seqmap) = synthetic_seqmap(&dir);

    // The edit: remove one A from the AAAAAA run at interbase [3, 9).
    //   MAF form (trimmed, `-` placeholder): delete the A at 1-based 5.
    //   VCF form (anchored):                 POS=4 REF=AA ALT=A.
    let del_maf = dir.join("del.maf");
    std::fs::write(
        &del_maf,
        format!(
            "{MAF_HEADER}\n\
             SYN\t0\tSage\tSYN1\ttinychr\t5\t5\t+\tframeshift_variant\tFrame_Shift_Del\tDEL\t\
             A\tA\t-\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t10\t20\n"
        ),
    )
    .unwrap();

    let (maf_alleles, maf_obs, stderr) = run_maf(
        &dir,
        &del_maf,
        &seqmap,
        &["--reference", fasta.to_str().unwrap()],
    );
    assert!(
        !stderr.contains("NOT fully justified"),
        "a reference was supplied, so nothing should be flagged: {stderr}"
    );
    let maf_allele = &json_lines(&maf_alleles)[0];
    assert_eq!(
        id_of(maf_allele),
        GOLDEN_DEL_ID,
        "MAF deletion must reach the vrs-python cross-validated id, got {maf_allele}"
    );
    // Fully-justified indels are RLEs, exactly as on the VCF path.
    assert_eq!(maf_allele["state"]["type"], "ReferenceLengthExpression");
    assert_eq!(maf_allele["state"]["length"], 5);
    assert_eq!(maf_allele["state"]["repeatSubunitLength"], 1);
    // The MAF-form alleles survive on the observation for debugging/round-tripping.
    let o = &json_lines(&maf_obs)[0];
    assert_eq!(o["referenceBases"], "A");
    assert_eq!(o["alternateBases"], "-");
    assert_eq!(o["sourcePos"], 5);
    assert_eq!(o["variantType"], "DEL");

    // Same edit through the VCF front end → same id.
    let vcf = dir.join("del.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         tinychr\t4\t.\tAA\tA\t.\tPASS\t.\n",
    )
    .unwrap();
    let vcf_alleles = dir.join("vcf_alleles.ndjson");
    let status = bin()
        .args(["--vcf"])
        .arg(&vcf)
        .arg("--seqmap")
        .arg(&seqmap)
        .arg("--reference")
        .arg(&fasta)
        .arg("--out-alleles")
        .arg(&vcf_alleles)
        .arg("--out-observations")
        .arg(dir.join("vcf_obs.ndjson"))
        .status()
        .unwrap();
    assert!(status.success());
    let vcf_allele = &json_lines(&std::fs::read_to_string(&vcf_alleles).unwrap())[0];
    assert_eq!(
        id_of(maf_allele),
        id_of(vcf_allele),
        "MAF and VCF front ends disagree on the same deletion"
    );

    // And the mirror case: an insertion of one A into the same run.
    //   MAF form: ref `-`, alt A, Start=5, End=6 (insertion sits *after* base 5).
    //   VCF form: POS=4 REF=A ALT=AA.
    let ins_maf = dir.join("ins.maf");
    std::fs::write(
        &ins_maf,
        format!(
            "{MAF_HEADER}\n\
             SYN\t0\tSage\tSYN1\ttinychr\t5\t6\t+\tframeshift_variant\tFrame_Shift_Ins\tINS\t\
             -\t-\tA\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t10\t20\n"
        ),
    )
    .unwrap();
    let ins_dir = dir.join("ins");
    std::fs::create_dir_all(&ins_dir).unwrap();
    let (ins_alleles, _, _) = run_maf(
        &ins_dir,
        &ins_maf,
        &seqmap,
        &["--reference", fasta.to_str().unwrap()],
    );
    let ins_maf_allele = &json_lines(&ins_alleles)[0];

    let ins_vcf = dir.join("ins.vcf");
    std::fs::write(
        &ins_vcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         tinychr\t4\t.\tA\tAA\t.\tPASS\t.\n",
    )
    .unwrap();
    let ins_vcf_alleles = dir.join("ins_vcf_alleles.ndjson");
    let status = bin()
        .args(["--vcf"])
        .arg(&ins_vcf)
        .arg("--seqmap")
        .arg(&seqmap)
        .arg("--reference")
        .arg(&fasta)
        .arg("--out-alleles")
        .arg(&ins_vcf_alleles)
        .arg("--out-observations")
        .arg(dir.join("ins_vcf_obs.ndjson"))
        .status()
        .unwrap();
    assert!(status.success());
    let ins_vcf_allele = &json_lines(&std::fs::read_to_string(&ins_vcf_alleles).unwrap())[0];
    assert_eq!(
        id_of(ins_maf_allele),
        id_of(ins_vcf_allele),
        "MAF and VCF front ends disagree on the same insertion"
    );
    assert_ne!(id_of(ins_maf_allele), GOLDEN_DEL_ID);

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn indel_without_reference_is_flagged_not_silently_wrong() {
    let dir = scratch("noref");
    let (_fasta, seqmap) = synthetic_seqmap(&dir);
    let maf = dir.join("del.maf");
    std::fs::write(
        &maf,
        format!(
            "{MAF_HEADER}\n\
             SYN\t0\tSage\tSYN1\ttinychr\t5\t5\t+\tframeshift_variant\tFrame_Shift_Del\tDEL\t\
             A\tA\t-\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t10\t20\n\
             SYN\t0\tSage\tSYN1\ttinychr\t1\t1\t+\tmissense_variant\tMissense_Mutation\tSNP\t\
             C\tC\tT\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t10\t20\n"
        ),
    )
    .unwrap();

    let (alleles, _, stderr) = run_maf(&dir, &maf, &seqmap, &[]);
    let alleles = json_lines(&alleles);
    let del = alleles
        .iter()
        .find(|a| a["state"]["sequence"] == "")
        .expect("the deletion allele");
    // Without a reference the id is NOT the canonical one, and the record admits it.
    assert_eq!(del["fullyJustified"], false);
    assert_ne!(id_of(del), GOLDEN_DEL_ID);
    assert!(
        stderr.contains("1 indel allele(s) are NOT fully justified"),
        "the run summary must count un-justified indels: {stderr}"
    );
    // Substitutions are exact either way, so they carry no caveat flag.
    let snv = alleles
        .iter()
        .find(|a| a["state"]["sequence"] == "T")
        .expect("the SNV allele");
    assert!(snv.get("fullyJustified").is_none());

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn min_tumor_alt_count_filters_unsupported_calls() {
    let dir = scratch("minalt");
    let (_fasta, seqmap) = synthetic_seqmap(&dir);
    // 28% of `nfib_ctf_biobank_2025` rows have t_alt_count = 0 — no read in the tumor
    // supports the allele. Ingesting those as "this specimen carries this variant" is
    // wrong, but the filter must be opt-in so nothing is dropped silently.
    let maf = dir.join("depth.maf");
    std::fs::write(
        &maf,
        format!(
            "{MAF_HEADER}\n\
             SYN\t0\tSage\tSYN1\ttinychr\t1\t1\t+\tmissense_variant\tMissense_Mutation\tSNP\t\
             C\tC\tT\t\tSUPPORTED\tSAMPLE-N\t\t\t28\t8\t20\n\
             SYN\t0\tSage\tSYN1\ttinychr\t3\t3\t+\tmissense_variant\tMissense_Mutation\tSNP\t\
             T\tT\tG\t\tUNSUPPORTED\tSAMPLE-N\t\t\t28\t28\t0\n"
        ),
    )
    .unwrap();

    // Default: nothing is filtered.
    let (_, obs, stderr) = run_maf(&dir, &maf, &seqmap, &[]);
    assert_eq!(json_lines(&obs).len(), 2);
    assert!(!stderr.contains("t_alt_count <"), "{stderr}");

    let filtered = dir.join("filtered");
    std::fs::create_dir_all(&filtered).unwrap();
    let (alleles, obs, stderr) = run_maf(&filtered, &maf, &seqmap, &["--min-tumor-alt-count", "1"]);
    let obs = json_lines(&obs);
    assert_eq!(obs.len(), 1, "only the supported call survives");
    assert_eq!(obs[0]["tumorSampleBarcode"], "SUPPORTED");
    // The unsupported row's allele must not be minted either.
    assert_eq!(json_lines(&alleles).len(), 1);
    assert!(stderr.contains("1 rows dropped for t_alt_count < 1"), "{stderr}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn rows_without_a_vrs_identity_are_kept_as_unnormalized() {
    let dir = scratch("unnorm");
    let (_fasta, seqmap) = synthetic_seqmap(&dir);
    // Row 1: a contig absent from the seqmap. Row 2: the right contig on the WRONG
    // assembly — hashing it against a SYN1 refget accession would mint a wrong id.
    // Row 3: both tumor alleles equal the reference, i.e. no variant at all.
    let maf = dir.join("odd.maf");
    std::fs::write(
        &maf,
        format!(
            "{MAF_HEADER}\n\
             SYN\t0\tSage\tSYN1\tchrZ\t100\t100\t+\tmissense_variant\tMissense_Mutation\tSNP\t\
             C\tC\tT\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t10\t20\n\
             SYN\t0\tSage\tGRCh37\ttinychr\t1\t1\t+\tmissense_variant\tMissense_Mutation\tSNP\t\
             C\tC\tT\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t10\t20\n\
             SYN\t0\tSage\tSYN1\ttinychr\t1\t1\t+\tsynonymous_variant\tSilent\tSNP\t\
             C\tC\tC\t\tSAMPLE-A\tSAMPLE-N\t\t\t30\t30\t0\n"
        ),
    )
    .unwrap();

    let (alleles, obs, stderr) = run_maf(
        &dir,
        &maf,
        &seqmap,
        &[
            "--study-id",
            "test_study",
            "--source",
            "syn9/odd.maf",
            // A local id has no correct default namespace, so the caller names one.
            "--variant-id-prefix",
            "nf:variant/",
        ],
    );
    let alleles = json_lines(&alleles);
    assert_eq!(alleles.len(), 2, "both odd rows are kept, got {alleles:#?}");
    for a in &alleles {
        assert_eq!(a["type"], "UnnormalizedVariant");
        assert_eq!(a["unnormalized"], true);
    }
    // Deterministic local key from {assembly}:{chrom}:{pos}:{ref}:{alt} (issue #95), so
    // the same row from another study collapses to the same node.
    assert_eq!(alleles[0]["id"], "nf:variant/SYN1:chrZ:100:C:T");
    assert!(alleles[0]["reason"].as_str().unwrap().contains("seqmap"));
    assert_eq!(alleles[1]["id"], "nf:variant/GRCh37:tinychr:1:C:T");
    assert!(alleles[1]["reason"].as_str().unwrap().contains("GRCh37"));
    // Rejection costs the row its VRS id, not its provenance: each odd row still has an
    // observation naming the sample/study/source, pointed at the unnormalized node. The
    // reference-only row is not a variant at all, so it has none.
    let obs = json_lines(&obs);
    assert_eq!(obs.len(), 2, "sample provenance must survive rejection: {obs:#?}");
    for (o, a) in obs.iter().zip(&alleles) {
        assert_eq!(o["variant"], a["id"]);
        assert_eq!(o["biosample"], "SAMPLE-A");
        assert_eq!(o["matchedNormalSampleBarcode"], "SAMPLE-N");
        assert_eq!(o["studyId"], "test_study");
        assert_eq!(o["center"], "Sage");
        assert_eq!(o["sourceFile"], "syn9/odd.maf");
        assert_eq!(o["tumorAltCount"], 20);
    }
    // Unknown contig: no seqmap entry to take `referenceName` from, so the MAF's own
    // `Chromosome` stands in; the assembly mismatch row keeps the MAF's build.
    assert_eq!(obs[0]["sourceContig"], "chrZ");
    assert_eq!(obs[0]["referenceName"], "chrZ");
    assert_eq!(obs[0]["assemblyId"], "SYN1");
    assert_eq!(obs[1]["referenceName"], "tinychr");
    assert_eq!(obs[1]["assemblyId"], "GRCh37");
    assert!(
        stderr.contains("2 unnormalized variant(s) kept with 2 observation(s)"),
        "{stderr}"
    );
    assert!(stderr.contains("1 rows had no non-reference tumor allele"), "{stderr}");

    // --strict turns the same input into a hard failure instead.
    let out = bin()
        .args(["maf", "--maf"])
        .arg(&maf)
        .arg("--seqmap")
        .arg(&seqmap)
        .arg("--out-alleles")
        .arg(dir.join("strict_alleles.ndjson"))
        .arg("--out-observations")
        .arg(dir.join("strict_obs.ndjson"))
        .arg("--strict")
        .output()
        .unwrap();
    assert!(!out.status.success(), "--strict must fail on an odd row");

    let _ = std::fs::remove_dir_all(&dir);
}
