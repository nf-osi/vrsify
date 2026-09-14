//! VCF ingestion → VRS alleles + per-sample observations.
//!
//! I/O uses `noodles-vcf`, which transparently handles plain-text and bgzipped
//! (`.vcf.gz`) input, spec-conformant header/sample parsing, and record tokenization.
//! The per-record projection helpers here (`build_allele`, `zygosity_for`,
//! `is_concrete_alt`, annotation parsing) operate on the string values noodles yields.
//!
//! Scope / limitations:
//!   * Coordinates default to the *naive* VCF→VRS projection: `start = POS-1`
//!     (interbase), `end = POS-1+len(REF)`, `state = alt` — correct for SNVs and for
//!     already left-aligned input. Run `bcftools norm -m- -f REF` upstream. The
//!     reference-based "fully justified" VRS normalization (in `normalize.rs`) is
//!     applied when a reference FASTA is supplied.
//!   * Symbolic/structural ALTs (`<DEL>`, breakends) are skipped and counted.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{bail, ensure, Context, Result};

use crate::vrs::{Allele, Coordinate, SequenceLocation, SequenceReference, State};

/// Parse the VEP `CSQ` column layout out of the INFO field's *description* string, which
/// VEP writes as `Consequence annotations from Ensembl VEP. Format: Allele|Consequence|
/// IMPACT|SYMBOL|Gene|...`. Returns the ordered column names, or `None` if no `Format:`
/// marker is present. snpEff's `ANN` has a fixed layout and needs no header lookup.
pub fn csq_format_from_description(description: &str) -> Option<Vec<String>> {
    let idx = description.find("Format:")?;
    let rest = description[idx + "Format:".len()..].trim();
    Some(rest.split('|').map(|s| s.trim().to_string()).collect())
}

/// Per-contig reference info. `refget` is the VRS `SQ.<digest>` accession used for
/// variant identity; `assembly`/`reference_name` are the denormalized fields Beacon
/// clients query by.
#[derive(Debug, Clone)]
pub struct SeqInfo {
    pub refget: String,
    pub assembly: Option<String>,
    pub reference_name: String,
}

pub type SeqMap = HashMap<String, SeqInfo>;

/// Load a TSV seqmap. Columns: `contig <tab> refgetAccession [<tab> assemblyId
/// [<tab> referenceName]]`. Lines starting with `#` and blanks are ignored.
pub fn load_seqmap<R: BufRead>(reader: R) -> Result<SeqMap> {
    let mut map = SeqMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            bail!("seqmap line {}: need at least contig<tab>refgetAccession", i + 1);
        }
        let contig = cols[0].to_string();
        let reference_name = cols.get(3).map(|s| s.to_string()).unwrap_or_else(|| contig.clone());
        map.insert(
            contig,
            SeqInfo {
                refget: cols[1].to_string(),
                assembly: cols.get(2).filter(|s| !s.is_empty()).map(|s| s.to_string()),
                reference_name,
            },
        );
    }
    Ok(map)
}

/// Resolve a contig name against a seqmap, tolerating a `chr`-prefix mismatch between
/// the source file and the seqmap.
///
/// A seqmap derived from a UCSC-style GRCh38 FASTA is keyed `chr1`, `chr1_KI270706v1_random`,
/// while cBioPortal MAFs write `1`, `1_KI270706v1_random`. Returns the resolved seqmap
/// key alongside the entry, so callers can look the *same* name up in a FASTA.
pub fn resolve_contig<'a>(map: &'a SeqMap, name: &str) -> Option<(&'a str, &'a SeqInfo)> {
    if let Some((k, v)) = map.get_key_value(name) {
        return Some((k.as_str(), v));
    }
    if let Some((k, v)) = map.get_key_value(&format!("chr{name}")) {
        return Some((k.as_str(), v));
    }
    let stripped = name.strip_prefix("chr")?;
    map.get_key_value(stripped).map(|(k, v)| (k.as_str(), v))
}

/// Build a VRS Allele from a (contig, 1-based POS, REF, ALT) tuple, naive projection.
/// Correct for SNVs and already-left-aligned input; for indels in repeats use
/// [`build_allele_normalized`] with a reference.
pub fn build_allele(seq: &SeqInfo, pos_1based: i64, reference: &str, alt: &str) -> Allele {
    let start = pos_1based - 1;
    let end = start + reference.len() as i64;
    allele_from_interval(seq, start, end, alt)
}

fn allele_from_interval(seq: &SeqInfo, start: i64, end: i64, alt: &str) -> Allele {
    Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: seq.refget.clone(),
            },
            start: Some(Coordinate::Definite(start)),
            end: Some(Coordinate::Definite(end)),
        },
        state: State::LiteralSequenceExpression {
            sequence: alt.to_string(),
        },
    }
}

/// Build a VRS Allele with reference-based fully-justified normalization. Applies
/// [`crate::normalize::normalize`] against `reference` so that equivalent indel
/// representations in repeat regions collapse to one canonical VRS id, matching
/// vrs-python. Falls back to the naive projection semantics for SNVs (unchanged).
/// Returns an error if the interval, REF bases, or seqmap accession disagree with
/// the supplied reference.
pub fn build_allele_normalized(
    seq: &SeqInfo,
    pos_1based: i64,
    reference: &str,
    alt: &str,
    ref_seq: &dyn crate::normalize::Reference,
) -> Result<Allele> {
    use crate::normalize::NormalizedState;
    ensure!(pos_1based >= 1, "POS must be >= 1");
    let start = usize::try_from(pos_1based - 1).context("POS is too large")?;
    let end = start
        .checked_add(reference.len())
        .context("REF interval overflow")?;
    crate::normalize::validate_reference(ref_seq, start, end, reference.as_bytes(), &seq.refget)?;
    let n = crate::normalize::normalize(ref_seq, start, end, alt.as_bytes())?;
    let alt_str = std::str::from_utf8(&n.alt).unwrap_or("").to_string();
    let state = match n.state {
        NormalizedState::Literal => State::LiteralSequenceExpression { sequence: alt_str },
        NormalizedState::ReferenceLengthExpression {
            length,
            repeat_subunit_length,
        } => State::ReferenceLengthExpression {
            length: length as i64,
            repeat_subunit_length: repeat_subunit_length as i64,
            // vrs-python default rle_seq_limit=50: include the sequence only when short.
            sequence: (alt_str.len() <= 50).then_some(alt_str),
        },
    };
    Ok(Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: seq.refget.clone(),
            },
            start: Some(Coordinate::Definite(n.start as i64)),
            end: Some(Coordinate::Definite(n.end as i64)),
        },
        state,
    })
}

/// True if an ALT string is a concrete sequence (not symbolic / breakend / missing).
pub fn is_concrete_alt(alt: &str) -> bool {
    !(alt.is_empty()
        || alt == "*"
        || alt == "."
        || alt.starts_with('<')
        || alt.contains('[')
        || alt.contains(']'))
}

/// Functional annotation extracted from a VCF INFO field (VEP `CSQ` / snpEff `ANN`).
/// Backs Beacon `geneId` / `aminoacidChange` queries. Context-full, so it lives on
/// the observation, never on the (context-free) VRS allele.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Annotation {
    /// Gene symbol (VEP `SYMBOL`, snpEff gene name field).
    pub gene_symbol: Option<String>,
    /// Stable gene id (VEP `Gene` = Ensembl id, snpEff gene id field).
    pub gene_id: Option<String>,
    /// Protein/amino-acid change (VEP `HGVSp` / `Amino_acids`, snpEff HGVS.p field).
    pub aminoacid_change: Option<String>,
    /// Molecular consequence term(s) (e.g. `missense_variant`).
    pub consequence: Option<String>,
}

impl Annotation {
    pub fn is_empty(&self) -> bool {
        self.gene_symbol.is_none()
            && self.gene_id.is_none()
            && self.aminoacid_change.is_none()
            && self.consequence.is_none()
    }
}

/// Parse a raw VCF INFO string (`k=v;flag;k2=v2`) into a key→value map. Flags map to
/// an empty string. Robust to spec-conformant INFO fields.
pub fn parse_info(info: &str) -> HashMap<&str, &str> {
    let mut m = HashMap::new();
    if info.is_empty() || info == "." {
        return m;
    }
    for entry in info.split(';') {
        if entry.is_empty() {
            continue;
        }
        match entry.split_once('=') {
            Some((k, v)) => {
                m.insert(k, v);
            }
            None => {
                m.insert(entry, "");
            }
        }
    }
    m
}

/// Extract a gene/consequence annotation from a raw INFO string, trying VEP `CSQ`
/// first, then snpEff `ANN`. Only the first (most-severe, per convention) transcript
/// annotation is used. Returns an empty `Annotation` when neither field is present.
///
/// VEP `CSQ` format string lives in the header `##INFO=<ID=CSQ,...Format: A|B|C>`; when
/// `csq_format` is supplied we key columns by name (`Consequence`, `SYMBOL`, `Gene`,
/// `HGVSp`/`Amino_acids`). Without it we fall back to VEP's default column order.
/// snpEff `ANN` has a fixed column order (Sequence Ontology `ANN` spec).
pub fn extract_annotation(info: &str, csq_format: Option<&[String]>) -> Annotation {
    let map = parse_info(info);
    if let Some(csq) = map.get("CSQ").or_else(|| map.get("vep")) {
        if let Some(first) = csq.split(',').next() {
            return parse_vep_csq(first, csq_format);
        }
    }
    if let Some(ann) = map.get("ANN") {
        if let Some(first) = ann.split(',').next() {
            return parse_snpeff_ann(first);
        }
    }
    Annotation::default()
}

/// Default VEP CSQ column order (the common subset), used when the header Format string
/// is unavailable. VEP's true order is header-defined; supply `csq_format` for fidelity.
const VEP_DEFAULT_FIELDS: &[&str] = &[
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "HGVSc",
    "HGVSp",
];

fn parse_vep_csq(entry: &str, csq_format: Option<&[String]>) -> Annotation {
    let cols: Vec<&str> = entry.split('|').collect();
    let col = |name: &str| -> Option<String> {
        let idx = match csq_format {
            Some(fmt) => fmt.iter().position(|f| f == name),
            None => VEP_DEFAULT_FIELDS.iter().position(|f| *f == name),
        }?;
        cols.get(idx)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(String::from)
    };
    let aa = col("HGVSp").or_else(|| col("Amino_acids"));
    Annotation {
        gene_symbol: col("SYMBOL"),
        gene_id: col("Gene"),
        aminoacid_change: aa,
        consequence: col("Consequence"),
    }
}

// snpEff ANN column order (fixed by the ANN spec):
// Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type |
// Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | ...
fn parse_snpeff_ann(entry: &str) -> Annotation {
    let c: Vec<&str> = entry.split('|').collect();
    let get = |i: usize| -> Option<String> {
        c.get(i)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(String::from)
    };
    Annotation {
        consequence: get(1),
        gene_symbol: get(3),
        gene_id: get(4),
        aminoacid_change: get(10),
    }
}

/// A per-sample observation of an allele.
#[derive(Debug, Clone)]
pub struct Observation {
    pub variant_id: String,
    pub sample: String,
    pub zygosity: String,
    pub contig: String,
    pub pos: i64,
    pub reference: String,
    pub alt: String,
    pub assembly: Option<String>,
    pub reference_name: String,
    pub source: String,
    pub annotation: Annotation,
}

impl Observation {
    pub fn to_json(&self) -> serde_json::Value {
        let mut m = serde_json::Map::new();
        m.insert("type".into(), "VariantCall".into());
        m.insert("variant".into(), self.variant_id.clone().into());
        m.insert("biosample".into(), self.sample.clone().into());
        m.insert("zygosity".into(), self.zygosity.clone().into());
        m.insert("referenceName".into(), self.reference_name.clone().into());
        if let Some(a) = &self.assembly {
            m.insert("assemblyId".into(), a.clone().into());
        }
        m.insert("sourceContig".into(), self.contig.clone().into());
        m.insert("sourcePos".into(), self.pos.into());
        m.insert("referenceBases".into(), self.reference.clone().into());
        m.insert("alternateBases".into(), self.alt.clone().into());
        m.insert("sourceFile".into(), self.source.clone().into());
        // Functional annotation (context-full → observation, not the allele).
        if let Some(g) = &self.annotation.gene_id {
            m.insert("affectedGene".into(), g.clone().into());
        }
        if let Some(s) = &self.annotation.gene_symbol {
            m.insert("affectedGeneSymbol".into(), s.clone().into());
        }
        if let Some(a) = &self.annotation.aminoacid_change {
            m.insert("aminoacidChange".into(), a.clone().into());
        }
        if let Some(c) = &self.annotation.consequence {
            m.insert("molecularConsequence".into(), c.clone().into());
        }
        serde_json::Value::Object(m)
    }
}

/// Given a sample GT string and the 1-based ALT allele number, return the zygosity if
/// the sample carries that ALT, else `None`. Handles phased (`|`) and unphased (`/`).
pub fn zygosity_for(gt: &str, alt_number: usize) -> Option<String> {
    let gt = gt.split(':').next().unwrap_or(gt); // GT is the first FORMAT subfield
    if gt == "." || gt.is_empty() {
        return None;
    }
    let alleles: Vec<Option<usize>> = gt
        .split(|c| c == '/' || c == '|')
        .map(|a| a.parse::<usize>().ok())
        .collect();
    let carries = alleles.iter().flatten().filter(|&&a| a == alt_number).count();
    if carries == 0 {
        return None;
    }
    let called = alleles.iter().flatten().count();
    Some(match (called, carries) {
        (1, _) => "hemizygous".to_string(),
        (n, c) if c == n => "homozygous".to_string(),
        _ => "heterozygous".to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zygosity() {
        assert_eq!(zygosity_for("0/1", 1).as_deref(), Some("heterozygous"));
        assert_eq!(zygosity_for("1|1", 1).as_deref(), Some("homozygous"));
        assert_eq!(zygosity_for("0/0", 1), None);
        assert_eq!(zygosity_for("1", 1).as_deref(), Some("hemizygous"));
        assert_eq!(zygosity_for("1/2", 2).as_deref(), Some("heterozygous"));
        assert_eq!(zygosity_for("./.", 1), None);
        assert_eq!(zygosity_for("0/1:35,40:75", 1).as_deref(), Some("heterozygous"));
    }

    #[test]
    fn vep_csq_with_header_format() {
        // VEP CSQ with an explicit header Format order.
        let fmt: Vec<String> = "Allele|Consequence|IMPACT|SYMBOL|Gene|HGVSp"
            .split('|')
            .map(String::from)
            .collect();
        let info = "AC=1;CSQ=T|missense_variant|MODERATE|APOE|ENSG00000130203|ENSP00000252486.3:p.Cys130Arg";
        let ann = extract_annotation(info, Some(&fmt));
        assert_eq!(ann.gene_symbol.as_deref(), Some("APOE"));
        assert_eq!(ann.gene_id.as_deref(), Some("ENSG00000130203"));
        assert_eq!(ann.consequence.as_deref(), Some("missense_variant"));
        assert_eq!(
            ann.aminoacid_change.as_deref(),
            Some("ENSP00000252486.3:p.Cys130Arg")
        );
    }

    #[test]
    fn vep_csq_default_order() {
        // No header format supplied → fall back to VEP default column order.
        let info = "CSQ=A|stop_gained|HIGH|BRCA1|ENSG00000012048|Transcript|ENST0|protein_coding|c.1|p.Arg100Ter";
        let ann = extract_annotation(info, None);
        assert_eq!(ann.gene_symbol.as_deref(), Some("BRCA1"));
        assert_eq!(ann.gene_id.as_deref(), Some("ENSG00000012048"));
        assert_eq!(ann.consequence.as_deref(), Some("stop_gained"));
        assert_eq!(ann.aminoacid_change.as_deref(), Some("p.Arg100Ter"));
    }

    #[test]
    fn snpeff_ann() {
        let info = "ANN=T|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST1|protein_coding|5/11|c.215C>G|p.Pro72Arg|215/1182|215/1182|72/393||";
        let ann = extract_annotation(info, None);
        assert_eq!(ann.gene_symbol.as_deref(), Some("TP53"));
        assert_eq!(ann.gene_id.as_deref(), Some("ENSG00000141510"));
        assert_eq!(ann.consequence.as_deref(), Some("missense_variant"));
        assert_eq!(ann.aminoacid_change.as_deref(), Some("p.Pro72Arg"));
    }

    #[test]
    fn no_annotation() {
        assert!(extract_annotation("AC=1;DP=30", None).is_empty());
        assert!(extract_annotation(".", None).is_empty());
    }

    #[test]
    fn resolve_contig_bridges_the_chr_prefix() {
        let mut map = SeqMap::new();
        let info = SeqInfo {
            refget: "SQ.test".into(),
            assembly: Some("GRCh38".into()),
            reference_name: "1".into(),
        };
        map.insert("chr1".into(), info.clone());
        // A MAF writes `1` where a UCSC-derived seqmap is keyed `chr1`; the resolved key
        // comes back so the caller can find the same contig in the FASTA.
        assert_eq!(resolve_contig(&map, "1").map(|(k, _)| k), Some("chr1"));
        assert_eq!(resolve_contig(&map, "chr1").map(|(k, _)| k), Some("chr1"));
        assert!(resolve_contig(&map, "2").is_none());
        // ...and the other direction, for an Ensembl-derived seqmap.
        let mut ens = SeqMap::new();
        ens.insert("17".into(), info);
        assert_eq!(resolve_contig(&ens, "chr17").map(|(k, _)| k), Some("17"));
        // Alt/random contigs keep their suffix through the prefix swap.
        let mut alt = SeqMap::new();
        alt.insert(
            "chr1_KI270706v1_random".into(),
            SeqInfo {
                refget: "SQ.alt".into(),
                assembly: None,
                reference_name: "1_KI270706v1_random".into(),
            },
        );
        assert_eq!(
            resolve_contig(&alt, "1_KI270706v1_random").map(|(k, _)| k),
            Some("chr1_KI270706v1_random")
        );
    }

    #[test]
    fn snv_allele_shape() {
        let seq = SeqInfo {
            refget: "SQ.test".into(),
            assembly: Some("GRCh38".into()),
            reference_name: "17".into(),
        };
        // A C>T SNV at 1-based POS 100 → interbase [99,100), state "T".
        let a = build_allele(&seq, 100, "C", "T");
        // ga4gh_serialize should reference location by digest and inline state.
        assert!(a.ga4gh_id().starts_with("ga4gh:VA."));
    }
}
