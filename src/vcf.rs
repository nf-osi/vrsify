//! VCF ingestion → VRS alleles + per-sample observations.
//!
//! I/O uses `noodles-vcf`, which transparently handles plain-text and bgzipped
//! (`.vcf.gz`) input, spec-conformant header/sample parsing, and record tokenization.
//! The per-record projection helpers here (`build_allele`, `zygosity_for`,
//! `is_concrete_alt`, annotation parsing) operate on the string values noodles yields.
//!
//! Scope / limitations:
//!   * Without a reference FASTA, the VCF→VRS projection trims the common REF/ALT
//!     suffix/prefix (the reference-free part of VRS normalization), so substitutions —
//!     including padded ones like `GT>GA` — are exact; pure indels still need the
//!     reference-based "fully justified" normalization (in `normalize.rs`), applied when
//!     a reference FASTA is supplied. Run `bcftools norm -m- -f REF` upstream as usual.
//!   * Symbolic/structural ALTs (`<DEL>`, breakends, `*`) are an explicit scope
//!     exclusion: they are counted and skipped, because their meaning lives in
//!     `INFO/END`/`SVLEN` rather than in the REF/ALT strings this module projects.
//!   * Non-ACGTN REF/ALT strings and contigs absent from the seqmap get no VRS id, but
//!     are *kept* as [`crate::unnormalized::UnnormalizedVariant`] nodes with their
//!     observations (in `main.rs`, mirroring the MAF front end).

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

/// Build a VRS Allele from a (contig, 1-based POS, REF, ALT) tuple without a reference
/// sequence, applying the reference-free part of VRS normalization: common suffix/prefix
/// trimming ([`crate::normalize::trim_common`]) plus the identity-allele RLE state. The
/// returned bool says whether the id is exact ("fully justified"): true for
/// substitutions and identity alleles — a padded `GT>GA` reaches the same id as the
/// reference-based path — false for pure indels, which need a reference to justify
/// across repeats (use [`build_allele_normalized`]).
pub fn build_allele(seq: &SeqInfo, pos_1based: i64, reference: &str, alt: &str) -> (Allele, bool) {
    // VCF nucleotide strings are case-insensitive (VCF 4.x sec. 1.6.1: bases may be in
    // either case), so `A>t` and `A>T` are the same edit and must hash to one id. The
    // reference-based path gets this from `normalize`, which uppercases internally.
    let reference = reference.to_ascii_uppercase();
    let alt = alt.to_ascii_uppercase();
    let start = pos_1based - 1;
    let end = start + reference.len() as i64;

    // Identity allele: the same RLE state `normalize` chooses, so the id agrees with the
    // reference-based path (which is exact here — no indel unit moves).
    if reference == alt {
        let state = State::ReferenceLengthExpression {
            length: alt.len() as i64,
            repeat_subunit_length: alt.len() as i64,
            sequence: (alt.len() <= 50).then_some(alt),
        };
        return (allele_at(seq, start, end, state), true);
    }

    let (pfx, sfx) = crate::normalize::trim_common(reference.as_bytes(), alt.as_bytes());
    let trimmed_alt = String::from_utf8_lossy(&alt.as_bytes()[pfx..alt.len() - sfx]).to_string();
    // Both sides non-empty after trimming → substitution/MNV, which VRS does not roll:
    // the trimmed literal IS the canonical form. One side empty → pure indel, not exact.
    let is_substitution = reference.len() > pfx + sfx && !trimmed_alt.is_empty();
    let allele = allele_at(
        seq,
        start + pfx as i64,
        end - sfx as i64,
        State::LiteralSequenceExpression { sequence: trimmed_alt },
    );
    (allele, is_substitution)
}

fn allele_at(seq: &SeqInfo, start: i64, end: i64, state: State) -> Allele {
    Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: seq.refget.clone(),
            },
            start: Some(Coordinate::Definite(start)),
            end: Some(Coordinate::Definite(end)),
        },
        state,
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

/// Extract a gene/consequence annotation for **one ALT allele** from a raw INFO string,
/// trying VEP `CSQ` first, then snpEff `ANN`. Returns an empty `Annotation` when neither
/// field is present.
///
/// `reference`/`alt` are the record's REF and the single ALT being observed, and
/// `alt_number` is that ALT's 1-based index in the record. A `CSQ`/`ANN` value holds one
/// entry per (allele, transcript) pair, so on a multiallelic record — or on a record
/// split by `bcftools norm -m-`, which leaves the *whole* original `CSQ` on every split
/// line — blindly taking the first entry can describe a different ALT than the one being
/// observed. Entries are therefore filtered to this ALT first (see
/// [`csq_entries_for_alt`]); the first survivor is used, which is VEP's most-severe
/// ordering. When no entry can be attributed to this ALT the whole list is used, so a
/// layout without an allele column degrades to the previous behaviour instead of
/// dropping the annotation.
///
/// VEP `CSQ` format string lives in the header `##INFO=<ID=CSQ,...Format: A|B|C>`; when
/// `csq_format` is supplied we key columns by name (`Consequence`, `SYMBOL`, `Gene`,
/// `HGVSp`/`Amino_acids`). Without it we fall back to VEP's default column order.
/// snpEff `ANN` has a fixed column order (Sequence Ontology `ANN` spec).
pub fn extract_annotation(
    info: &str,
    csq_format: Option<&[String]>,
    reference: &str,
    alt: &str,
    alt_number: usize,
) -> Annotation {
    let map = parse_info(info);
    if let Some(csq) = map.get("CSQ").or_else(|| map.get("vep")) {
        let allele_col = match csq_format {
            Some(fmt) => fmt.iter().position(|f| f == "Allele"),
            None => VEP_DEFAULT_FIELDS.iter().position(|f| *f == "Allele"),
        };
        let allele_num_col = csq_format.and_then(|f| f.iter().position(|c| c == "ALLELE_NUM"));
        let entries =
            csq_entries_for_alt(csq, allele_col, allele_num_col, reference, alt, alt_number);
        if let Some(first) = entries.first() {
            return parse_vep_csq(first, csq_format);
        }
    }
    if let Some(ann) = map.get("ANN") {
        // snpEff's ANN layout is fixed: column 0 is the allele, and there is no
        // ALLELE_NUM equivalent.
        let entries = csq_entries_for_alt(ann, Some(0), None, reference, alt, alt_number);
        if let Some(first) = entries.first() {
            return parse_snpeff_ann(first);
        }
    }
    Annotation::default()
}

/// The `CSQ`/`ANN` entries that describe `alt`, in source order.
///
/// Matching prefers an explicit `ALLELE_NUM` column (VEP `--allele_number`, the only
/// unambiguous signal) and otherwise compares the entry's allele column against the
/// spellings a tool may use for this ALT (see [`vep_allele_spellings`]). If nothing
/// matches — no allele column, an unrecognized spelling — every entry is returned, so
/// the caller falls back to the first-entry convention rather than losing the annotation.
fn csq_entries_for_alt<'a>(
    value: &'a str,
    allele_col: Option<usize>,
    allele_num_col: Option<usize>,
    reference: &str,
    alt: &str,
    alt_number: usize,
) -> Vec<&'a str> {
    let all: Vec<&str> = value.split(',').filter(|e| !e.is_empty()).collect();
    let field = |entry: &'a str, col: usize| entry.split('|').nth(col).map(str::trim);

    if let Some(col) = allele_num_col {
        let matched: Vec<&str> = all
            .iter()
            .copied()
            .filter(|e| field(e, col) == Some(alt_number.to_string().as_str()))
            .collect();
        if !matched.is_empty() {
            return matched;
        }
    }
    if let Some(col) = allele_col {
        let spellings = vep_allele_spellings(reference, alt);
        let matched: Vec<&str> = all
            .iter()
            .copied()
            .filter(|e| {
                field(e, col)
                    .is_some_and(|a| spellings.iter().any(|s| a.eq_ignore_ascii_case(s)))
            })
            .collect();
        if !matched.is_empty() {
            return matched;
        }
    }
    all
}

/// The spellings a VEP/snpEff `Allele` column may use for one VCF REF/ALT pair.
///
/// snpEff writes the ALT verbatim. VEP writes the *minimal* allele: for an indel it
/// strips the bases REF and ALT share on the left and writes `-` when nothing is left,
/// so a VCF `CTT>CT` deletion is reported as `T` — sometimes as `-`, depending on how
/// much padding the caller used. All of these are accepted.
fn vep_allele_spellings(reference: &str, alt: &str) -> Vec<String> {
    let mut out = vec![alt.to_string()];
    if reference.len() != alt.len() {
        // VEP trims the shared left flank only.
        let pfx = reference
            .bytes()
            .zip(alt.bytes())
            .take_while(|(r, a)| r.eq_ignore_ascii_case(a))
            .count();
        let trimmed = &alt[pfx..];
        out.push(if trimmed.is_empty() { "-".to_string() } else { trimmed.to_string() });
        // Some writers trim both flanks (`bcftools`-style minimal representation).
        let (pfx, sfx) = crate::normalize::trim_common(reference.as_bytes(), alt.as_bytes());
        let trimmed = &alt[pfx..alt.len() - sfx];
        out.push(if trimmed.is_empty() { "-".to_string() } else { trimmed.to_string() });
    }
    out
}

/// Decode the percent-escapes VCF 4.3 (sec. 1.2) requires in INFO values, which VEP uses
/// for characters that would otherwise break the field: `%3D` for `=` (an HGVSp
/// synonymous change is written `p.Cys130%3D`), `%2C` for `,`, `%3B` for `;`, `%25` for
/// `%`. Decoding happens *after* the value is split on `,` and `|`, so an escaped
/// delimiter cannot be mistaken for a real one. Malformed escapes are left as written.
fn percent_decode(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%'
            && i + 2 < bytes.len()
            && let (Some(hi), Some(lo)) = (
                (bytes[i + 1] as char).to_digit(16),
                (bytes[i + 2] as char).to_digit(16),
            )
        {
            out.push((hi * 16 + lo) as u8);
            i += 3;
            continue;
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
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
            .map(percent_decode)
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
            .map(percent_decode)
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
        // Issue #95's class for "sample S carries variant V". Both front ends label
        // the relationship the same way, so one loader handles either stream.
        m.insert("type".into(), "VariantObservation".into());
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
///
/// Ploidy comes from the number of GT fields, not from the number of *called* ones: a
/// partially missing call (`1/.`, VCF 4.x sec. 1.6.2) is a diploid genotype whose second
/// allele is unknown, so it could be `1/0` or `1/1`. Reporting `hemizygous` there would
/// assert a single-copy locus (chrX/chrY in a male sample, a haploid `1` call) that the
/// data does not support, so the uncertainty is kept as `"unknown"`.
pub fn zygosity_for(gt: &str, alt_number: usize) -> Option<String> {
    let gt = gt.split(':').next().unwrap_or(gt); // GT is the first FORMAT subfield
    if gt == "." || gt.is_empty() {
        return None;
    }
    let alleles: Vec<Option<usize>> = gt
        .split(['/', '|'])
        .map(|a| a.parse::<usize>().ok())
        .collect();
    let carries = alleles.iter().flatten().filter(|&&a| a == alt_number).count();
    if carries == 0 {
        return None;
    }
    let ploidy = alleles.len();
    let called = alleles.iter().flatten().count();
    if called < ploidy {
        // The sample carries the ALT, but the copy count is not knowable from this GT.
        return Some("unknown".to_string());
    }
    Some(match (ploidy, carries) {
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
        // Partially missing: diploid with one unknown allele, NOT hemizygous.
        assert_eq!(zygosity_for("1/.", 1).as_deref(), Some("unknown"));
        assert_eq!(zygosity_for("./1", 1).as_deref(), Some("unknown"));
        assert_eq!(zygosity_for(".|1", 1).as_deref(), Some("unknown"));
        assert_eq!(zygosity_for("1/1/.", 1).as_deref(), Some("unknown"));
        assert_eq!(zygosity_for("0/.", 1), None);
        assert_eq!(zygosity_for("1/1/1", 1).as_deref(), Some("homozygous"));
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
        let ann = extract_annotation(info, Some(&fmt), "C", "T", 1);
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
        let ann = extract_annotation(info, None, "C", "A", 1);
        assert_eq!(ann.gene_symbol.as_deref(), Some("BRCA1"));
        assert_eq!(ann.gene_id.as_deref(), Some("ENSG00000012048"));
        assert_eq!(ann.consequence.as_deref(), Some("stop_gained"));
        assert_eq!(ann.aminoacid_change.as_deref(), Some("p.Arg100Ter"));
    }

    #[test]
    fn snpeff_ann() {
        let info = "ANN=T|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST1|protein_coding|5/11|c.215C>G|p.Pro72Arg|215/1182|215/1182|72/393||";
        let ann = extract_annotation(info, None, "C", "T", 1);
        assert_eq!(ann.gene_symbol.as_deref(), Some("TP53"));
        assert_eq!(ann.gene_id.as_deref(), Some("ENSG00000141510"));
        assert_eq!(ann.consequence.as_deref(), Some("missense_variant"));
        assert_eq!(ann.aminoacid_change.as_deref(), Some("p.Pro72Arg"));
    }

    #[test]
    fn no_annotation() {
        assert!(extract_annotation("AC=1;DP=30", None, "C", "T", 1).is_empty());
        assert!(extract_annotation(".", None, "C", "T", 1).is_empty());
    }

    /// A multiallelic record carries one CSQ entry per (allele, transcript). Taking the
    /// first entry unconditionally annotated ALT 2 with ALT 1's gene/consequence.
    #[test]
    fn csq_entry_is_matched_to_the_alt_being_observed() {
        let fmt: Vec<String> = "Allele|Consequence|IMPACT|SYMBOL|Gene|HGVSp"
            .split('|')
            .map(String::from)
            .collect();
        let info = "CSQ=T|missense_variant|MODERATE|APOE|ENSG1|p.Cys130Arg,\
                    G|stop_gained|HIGH|APOE|ENSG1|p.Cys130Ter";
        let info = &info.replace(' ', "");
        for (alt, consequence, aa) in [
            ("T", "missense_variant", "p.Cys130Arg"),
            ("G", "stop_gained", "p.Cys130Ter"),
        ] {
            let ann = extract_annotation(info, Some(&fmt), "C", alt, 1);
            assert_eq!(ann.consequence.as_deref(), Some(consequence), "ALT={alt}");
            assert_eq!(ann.aminoacid_change.as_deref(), Some(aa), "ALT={alt}");
        }
    }

    /// `bcftools norm -m-` splits a multiallelic record but leaves the *whole* original
    /// CSQ on every split line, so the allele column is the only way to tell them apart.
    /// The several transcript entries for the matching allele stay in VEP's most-severe
    /// order, so the first survivor is still the one to use.
    #[test]
    fn csq_keeps_most_severe_entry_among_the_matching_allele() {
        let fmt: Vec<String> = "Allele|Consequence|IMPACT|SYMBOL|Gene"
            .split('|')
            .map(String::from)
            .collect();
        let info = "CSQ=T|missense_variant|MODERATE|APOE|ENSG1,\
                    T|intron_variant|MODIFIER|APOE|ENSG1,\
                    G|stop_gained|HIGH|OTHER|ENSG2";
        let info = &info.replace(' ', "");
        let ann = extract_annotation(info, Some(&fmt), "C", "T", 1);
        assert_eq!(ann.consequence.as_deref(), Some("missense_variant"));
        assert_eq!(ann.gene_symbol.as_deref(), Some("APOE"));
    }

    /// VEP's `ALLELE_NUM` names the ALT index outright; it wins over spelling matching.
    #[test]
    fn csq_prefers_allele_num_when_present() {
        let fmt: Vec<String> = "Allele|Consequence|SYMBOL|ALLELE_NUM"
            .split('|')
            .map(String::from)
            .collect();
        let info = "CSQ=T|missense_variant|APOE|1,T|stop_gained|OTHER|2";
        assert_eq!(
            extract_annotation(info, Some(&fmt), "C", "T", 2).gene_symbol.as_deref(),
            Some("OTHER")
        );
    }

    /// VEP reports the *minimal* allele, so a VCF-padded indel never matches the ALT
    /// verbatim: `CTT>CT` is `T` (or `-`) in the CSQ. Both spellings must resolve.
    #[test]
    fn csq_matches_veps_minimal_indel_allele() {
        let fmt: Vec<String> = "Allele|Consequence|SYMBOL".split('|').map(String::from).collect();
        for allele in ["T", "-"] {
            let info = format!("CSQ={allele}|frameshift_variant|NF1");
            let ann = extract_annotation(&info, Some(&fmt), "CTT", "CT", 1);
            assert_eq!(ann.gene_symbol.as_deref(), Some("NF1"), "CSQ Allele={allele}");
        }
    }

    /// With no attributable entry (a layout with no allele column, an unrecognized
    /// spelling) the annotation must degrade to the first entry, not vanish.
    #[test]
    fn csq_falls_back_to_the_first_entry_when_nothing_matches() {
        let fmt: Vec<String> = "Consequence|SYMBOL".split('|').map(String::from).collect();
        let info = "CSQ=missense_variant|APOE";
        assert_eq!(
            extract_annotation(info, Some(&fmt), "C", "T", 1).gene_symbol.as_deref(),
            Some("APOE")
        );
        // Allele column present but spelled in a way we do not recognize.
        let fmt: Vec<String> = "Allele|Consequence|SYMBOL".split('|').map(String::from).collect();
        let info = "CSQ=C/T|missense_variant|APOE";
        assert_eq!(
            extract_annotation(info, Some(&fmt), "C", "T", 1).gene_symbol.as_deref(),
            Some("APOE")
        );
    }

    /// VCF 4.3 percent-encodes characters that would break an INFO value; VEP writes a
    /// synonymous HGVSp as `p.Cys130%3D`, which must not reach consumers escaped.
    #[test]
    fn csq_fields_are_percent_decoded() {
        let fmt: Vec<String> = "Allele|Consequence|SYMBOL|HGVSp"
            .split('|')
            .map(String::from)
            .collect();
        let info = "CSQ=T|synonymous_variant|APOE|ENSP1:p.Cys130%3D";
        let ann = extract_annotation(info, Some(&fmt), "C", "T", 1);
        assert_eq!(ann.aminoacid_change.as_deref(), Some("ENSP1:p.Cys130="));
        assert_eq!(percent_decode("a%2Cb%3Bc%25d"), "a,b;c%d");
        // Malformed or truncated escapes are left exactly as written.
        assert_eq!(percent_decode("100%"), "100%");
        assert_eq!(percent_decode("%ZZ"), "%ZZ");
        assert_eq!(percent_decode("%3"), "%3");
        assert_eq!(percent_decode("plain"), "plain");
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

    fn test_seq() -> SeqInfo {
        SeqInfo {
            refget: "SQ.test".into(),
            assembly: Some("GRCh38".into()),
            reference_name: "17".into(),
        }
    }

    #[test]
    fn snv_allele_shape() {
        // A C>T SNV at 1-based POS 100 → interbase [99,100), state "T".
        let (a, exact) = build_allele(&test_seq(), 100, "C", "T");
        assert!(exact, "an SNV is exact without a reference");
        // ga4gh_serialize should reference location by digest and inline state.
        assert!(a.ga4gh_id().starts_with("ga4gh:VA."));
    }

    #[test]
    fn padded_substitution_trims_to_the_canonical_literal() {
        let seq = test_seq();
        // GT>GA at POS 100: the shared G pads the real T>A at 101. The trimmed form is
        // the canonical VRS literal, so the padded and minimal spellings share one id.
        let (padded, exact) = build_allele(&seq, 100, "GT", "GA");
        assert!(exact, "a substitution is exact after trimming");
        let (minimal, _) = build_allele(&seq, 101, "T", "A");
        assert_eq!(padded.ga4gh_id(), minimal.ga4gh_id());
        // Suffix pads trim too, shrinking `end`.
        let (suffix_padded, _) = build_allele(&seq, 101, "TG", "AG");
        assert_eq!(suffix_padded.ga4gh_id(), minimal.ga4gh_id());
    }

    #[test]
    fn identity_allele_matches_the_normalized_rle_state() {
        // REF == ALT: `normalize` emits an RLE spanning the original interval; the
        // no-reference path must mint the same id, not a LiteralSequenceExpression.
        let (a, exact) = build_allele(&test_seq(), 100, "AT", "AT");
        assert!(exact);
        let v = a.to_output_value();
        assert_eq!(v["state"]["type"], "ReferenceLengthExpression");
        assert_eq!(v["state"]["length"], 2);
        assert_eq!(v["state"]["repeatSubunitLength"], 2);
        assert_eq!(v["location"]["start"], 99);
        assert_eq!(v["location"]["end"], 101);
    }

    #[test]
    fn pure_indels_are_trimmed_but_not_exact() {
        let seq = test_seq();
        // Padded deletion AA>A: trims to a 1-base deletion, but justification across the
        // repeat needs a reference, so it is not exact.
        let (del, exact) = build_allele(&seq, 100, "AA", "A");
        assert!(!exact, "a pure indel cannot be exact without a reference");
        let v = del.to_output_value();
        assert_eq!(v["state"]["sequence"], "");
        assert_eq!(v["location"]["start"], 99);
        assert_eq!(v["location"]["end"], 100);
        let (_, exact) = build_allele(&seq, 100, "A", "AA");
        assert!(!exact);
    }
}
