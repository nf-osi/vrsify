//! Acceptance: the tool reads a *bgzipped* multi-sample VCF and produces the same
//! NDJSON it does for plain text (byte-exact rs7412 golden id, per-sample observations).

use std::io::Write;
use std::process::Command;

use noodles_bgzf as bgzf;

const VCF: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|HGVSp\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr19\t44908822\trs7412\tC\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|APOE|ENSG00000130203|p.Cys130Arg\tGT\t0/1\t1/1
chr19\t44908830\t.\tG\t<DEL>\t.\tPASS\t.\tGT\t0/1\t0/0
";

const SEQMAP: &str = "chr19\tSQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl\tGRCh38\t19\n";

#[test]
fn reads_bgzipped_multisample_vcf() {
    let dir = std::env::temp_dir().join(format!("vrsify_bgzf_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    // Write a bgzipped VCF (extension .vcf.gz so noodles selects the bgzf codec).
    let vcf_gz = dir.join("in.vcf.gz");
    {
        let f = std::fs::File::create(&vcf_gz).unwrap();
        let mut w = bgzf::io::Writer::new(f);
        w.write_all(VCF.as_bytes()).unwrap();
        w.finish().unwrap();
    }
    let seqmap = dir.join("seqmap.tsv");
    std::fs::write(&seqmap, SEQMAP).unwrap();

    let alleles = dir.join("alleles.ndjson");
    let obs = dir.join("obs.ndjson");

    let status = Command::new(env!("CARGO_BIN_EXE_vrsify"))
        .args([
            "--vcf",
            vcf_gz.to_str().unwrap(),
            "--seqmap",
            seqmap.to_str().unwrap(),
            "--out-alleles",
            alleles.to_str().unwrap(),
            "--out-observations",
            obs.to_str().unwrap(),
            "--source",
            "syn/test.vcf.gz",
        ])
        .status()
        .unwrap();
    assert!(status.success(), "vrsify exited non-zero");

    let alleles_txt = std::fs::read_to_string(&alleles).unwrap();
    assert!(
        alleles_txt.contains("ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"),
        "bgzipped input did not reproduce the rs7412 golden VA id:\n{alleles_txt}"
    );
    // Exactly one unique allele (the <DEL> ALT is skipped).
    assert_eq!(alleles_txt.lines().count(), 1);

    let obs_txt = std::fs::read_to_string(&obs).unwrap();
    // S1 heterozygous + S2 homozygous for the rs7412 ALT → two observations.
    assert_eq!(obs_txt.lines().count(), 2, "expected 2 observations\n{obs_txt}");
    assert!(obs_txt.contains("\"zygosity\":\"heterozygous\""));
    assert!(obs_txt.contains("\"zygosity\":\"homozygous\""));
    // VEP CSQ annotation carried onto observations.
    assert!(obs_txt.contains("\"affectedGeneSymbol\":\"APOE\""));
    assert!(obs_txt.contains("\"affectedGene\":\"ENSG00000130203\""));
    assert!(obs_txt.contains("\"aminoacidChange\":\"p.Cys130Arg\""));

    let _ = std::fs::remove_dir_all(&dir);
}
