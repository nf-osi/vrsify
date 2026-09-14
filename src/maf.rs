//! MAF (Mutation Annotation Format) ingestion → VRS alleles + somatic observations.
//!
//! The VRS identity engine (`digest`, `vrs`, `normalize`, `refget`) is format-agnostic;
//! only the *front end* is format-specific. This module is the MAF front end, the
//! sibling of [`crate::vcf`], and it exists because cBioPortal study data is
//! distributed as MAF (`data_mutations.txt`), not VCF.
//!
//! ## Why MAF → VRS is *simpler* than VCF → VRS
//!
//! VRS locations are **interbase** (0-based, half-open). MAF already stores indels in
//! trimmed form using `-` placeholders, so the projection is a pure coordinate shift
//! with **no anchor-base reference lookup** — unlike converting MAF to VCF, which must
//! fetch the flanking base to synthesise a VCF REF/ALT:
//!
//! | MAF row                                        | VRS interbase interval | state    |
//! |------------------------------------------------|------------------------|----------|
//! | SNP/DNP/ONP `Start=100 End=102 ref=ACG alt=TTT`| `[99, 102)`            | `"TTT"`  |
//! | DEL `Start=100 End=102 ref=ACG alt=-`          | `[99, 102)`            | `""`     |
//! | INS `Start=100 End=101 ref=- alt=TT`           | `[100, 100)`           | `"TT"`   |
//!
//! (An insertion sits *after* `Start_Position`; base `Start` occupies interbase
//! `[Start-1, Start)`, so the insertion point is the empty interval `[Start, Start)`.)
//!
//! ## Fully-justified normalization still needs a reference
//!
//! The table above is the *trimmed* projection. It is byte-exact for substitutions —
//! a MAF SNP and the same variant from a VCF yield the **same `ga4gh:VA.` id** — but
//! for indels in repeat regions VRS requires bidirectional ("fully justified")
//! justification against the reference, so `--reference` is required for indel ids to
//! agree with vrs-python, ClinVar, gnomAD, and the VCF path. Without it, indel alleles
//! are emitted with `"fullyJustified": false` and counted loudly.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{bail, ensure, Context, Result};
use serde_json::{Map, Value};

use crate::vcf::SeqInfo;
use crate::vrs::{Allele, Coordinate, SequenceLocation, SequenceReference, State};

/// Column names this module reads, with the spelling variants seen in the wild
/// (GDC/cBioPortal `Start_Position` vs. older TCGA `Start_position`).
const ALIASES: &[(&str, &[&str])] = &[
    ("Start_Position", &["Start_position", "start_position"]),
    ("End_Position", &["End_position", "end_position"]),
    ("Chromosome", &["chromosome", "Chrom", "Chr"]),
    ("NCBI_Build", &["ncbi_build", "Build"]),
];

/// A MAF header: column name → column index, tolerant of the known spelling variants.
#[derive(Debug, Clone)]
pub struct MafHeader {
    idx: HashMap<String, usize>,
}

/// Required columns. Everything else this module reads is optional.
const REQUIRED: &[&str] = &[
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
];

impl MafHeader {
    /// Parse a MAF header line (tab-separated column names).
    pub fn parse(line: &str) -> Result<Self> {
        let mut idx = HashMap::new();
        for (i, col) in line.trim_end_matches(['\r', '\n']).split('\t').enumerate() {
            idx.insert(col.trim().to_string(), i);
        }
        // Canonicalize known spelling variants onto the canonical name.
        for (canonical, variants) in ALIASES {
            if !idx.contains_key(*canonical)
                && let Some(i) = variants.iter().find_map(|v| idx.get(*v)).copied() {
                    idx.insert(canonical.to_string(), i);
                }
        }
        let header = Self { idx };
        let missing: Vec<&str> = REQUIRED
            .iter()
            .copied()
            .filter(|c| !header.idx.contains_key(*c))
            .collect();
        if !missing.is_empty() {
            bail!("MAF is missing required column(s): {}", missing.join(", "));
        }
        Ok(header)
    }

    /// Value of `column` in `row`, trimmed, with MAF's null spellings (`""`, `.`, `NA`)
    /// normalized to `None`. The `-` allele placeholder is *not* nulled here — callers
    /// need to distinguish it (see [`allele_bases`]).
    pub fn get<'a>(&self, row: &'a [&'a str], column: &str) -> Option<&'a str> {
        let v = row.get(*self.idx.get(column)?)?.trim();
        (!matches!(v, "" | "." | "NA")).then_some(v)
    }

    pub fn has(&self, column: &str) -> bool {
        self.idx.contains_key(column)
    }
}

/// Strip a MAF allele placeholder: `-` (and the null spellings) mean "no bases".
pub fn allele_bases(raw: Option<&str>) -> &str {
    match raw {
        Some("-") | None => "",
        Some(v) => v,
    }
}

/// The interbase interval + replacement sequence a MAF row projects to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MafInterval {
    /// Interbase (0-based, half-open) start.
    pub start: i64,
    /// Interbase end. Equals `start` for a pure insertion.
    pub end: i64,
    /// Replacement bases (uppercased). Empty for a pure deletion.
    pub alt: String,
    /// Reference bases the interval spans (uppercased). Empty for a pure insertion.
    pub reference: String,
}

impl MafInterval {
    pub fn is_indel(&self) -> bool {
        self.reference.len() != self.alt.len()
    }
}

/// Project a MAF row's coordinates/alleles into a VRS interbase interval.
///
/// `start_position`/`end_position` are MAF's 1-based inclusive coordinates;
/// `reference_allele`/`alt_allele` may be `-` placeholders. `end_position` is used only
/// as a consistency check — the `Reference_Allele` string is authoritative for the span,
/// because it is what the alleles are actually built from.
pub fn maf_to_interbase(
    start_position: i64,
    end_position: Option<i64>,
    reference_allele: &str,
    alt_allele: &str,
) -> Result<(MafInterval, Option<String>)> {
    // Accept either the raw MAF cell (`-`) or already-stripped bases.
    let reference = allele_bases(Some(reference_allele)).to_ascii_uppercase();
    let alt = allele_bases(Some(alt_allele)).to_ascii_uppercase();
    if start_position < 1 {
        bail!("Start_Position must be >= 1, got {start_position}");
    }

    let (interval, expected_end) = if reference.is_empty() {
        // Insertion: the variant sits *after* base `Start_Position`, i.e. in the empty
        // interbase interval [Start_Position, Start_Position). MAF convention sets
        // End_Position = Start_Position + 1 (the base on the other side).
        (
            MafInterval {
                start: start_position,
                end: start_position,
                alt,
                reference,
            },
            start_position
                .checked_add(1)
                .context("insertion coordinate overflow")?,
        )
    } else {
        // Substitution or deletion: `Reference_Allele` spans 1-based
        // [Start_Position, Start_Position + len - 1] → interbase [Start-1, Start-1+len).
        let len = reference.len() as i64;
        let end = (start_position - 1)
            .checked_add(len)
            .context("REF interval overflow")?;
        (
            MafInterval {
                start: start_position - 1,
                end,
                alt,
                reference,
            },
            end,
        )
    };

    // Report (don't fail on) a disagreeing End_Position: the allele string wins, but a
    // mismatch means the row is malformed and worth surfacing in the run summary.
    let warning = end_position.filter(|e| *e != expected_end).map(|e| {
        format!(
            "End_Position {e} disagrees with Start_Position {start_position} + \
             Reference_Allele '{reference_allele}' (expected {expected_end})"
        )
    });
    Ok((interval, warning))
}

/// The effective tumor alt allele for a MAF row.
///
/// `Tumor_Seq_Allele2` is the variant allele by convention; when it equals the reference
/// (some callers put the variant in allele1) fall back to `Tumor_Seq_Allele1`. Returns
/// `None` when neither differs from the reference, i.e. the row records no variant.
pub fn effective_alt<'a>(reference: &str, tsa2: &'a str, tsa1: &'a str) -> Option<&'a str> {
    if !tsa2.eq_ignore_ascii_case(reference) {
        return Some(tsa2);
    }
    if !tsa1.eq_ignore_ascii_case(reference) {
        return Some(tsa1);
    }
    None
}

/// True if `a` and `b` name the same assembly, tolerating the UCSC/GRC spellings that
/// MAF `NCBI_Build` and a FASTA-derived seqmap disagree on (`hg38` vs `GRCh38`).
pub fn build_matches(a: &str, b: &str) -> bool {
    canonical_build(a) == canonical_build(b)
}

fn canonical_build(build: &str) -> String {
    let b = build.trim().to_ascii_lowercase();
    match b.as_str() {
        "hg38" | "grch38" | "grch38.p13" | "38" => "grch38".to_string(),
        "hg19" | "grch37" | "37" => "grch37".to_string(),
        "hg18" | "ncbi36" | "36" => "ncbi36".to_string(),
        _ => b,
    }
}

/// Build a VRS allele for a MAF interval. With `reference` the allele is fully
/// justified; without it the trimmed projection is used, which is exact for
/// substitutions but not guaranteed canonical for indels in repeats. The returned bool
/// says whether the id is exact ("fully justified"): always with a reference, and for
/// substitutions without one — MAF alleles are trimmed by convention, but rows that
/// arrive padded are re-trimmed here so they reach the canonical id too.
/// With a reference, returns an error for interval, REF, or seqmap identity mismatches.
pub fn build_maf_allele(
    seq: &SeqInfo,
    interval: &MafInterval,
    reference: Option<&dyn crate::normalize::Reference>,
) -> Result<(Allele, bool)> {
    ensure!(
        interval.start >= 0 && interval.end >= interval.start,
        "invalid MAF interval"
    );
    Ok(match reference {
        Some(ref_seq) => {
            use crate::normalize::NormalizedState;
            let start = usize::try_from(interval.start).context("MAF start is too large")?;
            let end = usize::try_from(interval.end).context("MAF end is too large")?;
            crate::normalize::validate_reference(
                ref_seq,
                start,
                end,
                interval.reference.as_bytes(),
                &seq.refget,
            )?;
            let n = crate::normalize::normalize(
                ref_seq,
                start,
                end,
                interval.alt.as_bytes(),
            )?;
            let alt = String::from_utf8_lossy(&n.alt).to_string();
            let state = match n.state {
                NormalizedState::Literal => State::LiteralSequenceExpression { sequence: alt },
                NormalizedState::ReferenceLengthExpression {
                    length,
                    repeat_subunit_length,
                } => State::ReferenceLengthExpression {
                    length: length as i64,
                    repeat_subunit_length: repeat_subunit_length as i64,
                    // vrs-python default rle_seq_limit=50.
                    sequence: (alt.len() <= 50).then_some(alt),
                },
            };
            (allele_at(seq, n.start as i64, n.end as i64, state), true)
        }
        None => {
            ensure!(
                !(interval.reference.is_empty() && interval.alt.is_empty()),
                "cannot build an allele from an empty no-op edit"
            );
            // Identity allele: the same RLE state `normalize` chooses (exact — no indel
            // unit moves), so the id agrees with the reference-based path.
            if interval.reference == interval.alt {
                let state = State::ReferenceLengthExpression {
                    length: interval.alt.len() as i64,
                    repeat_subunit_length: interval.alt.len() as i64,
                    sequence: (interval.alt.len() <= 50).then_some(interval.alt.clone()),
                };
                return Ok((allele_at(seq, interval.start, interval.end, state), true));
            }
            // Reference-free trimming: exact for substitutions (the trimmed literal is
            // canonical); a pure indel still needs a reference to justify across repeats.
            let (pfx, sfx) = crate::normalize::trim_common(
                interval.reference.as_bytes(),
                interval.alt.as_bytes(),
            );
            let alt = String::from_utf8_lossy(
                &interval.alt.as_bytes()[pfx..interval.alt.len() - sfx],
            )
            .to_string();
            let is_substitution = interval.reference.len() > pfx + sfx && !alt.is_empty();
            let allele = allele_at(
                seq,
                interval.start + pfx as i64,
                interval.end - sfx as i64,
                State::LiteralSequenceExpression { sequence: alt },
            );
            (allele, is_substitution)
        }
    })
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

/// Functional annotation carried on a MAF row. Unlike the VCF path (which has to parse
/// VEP `CSQ` / snpEff `ANN` out of INFO) MAF has already flattened the picked
/// transcript into named columns, so this is a straight column read.
///
/// Context-full, so it belongs on the observation, never on the context-free allele.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MafAnnotation {
    pub gene_symbol: Option<String>,
    pub entrez_gene_id: Option<String>,
    pub hgnc_id: Option<String>,
    /// Ensembl gene id (MAF `Gene`).
    pub ensembl_gene_id: Option<String>,
    pub transcript_id: Option<String>,
    /// MAF `Variant_Classification` (MAF's own small vocabulary, e.g. `Missense_Mutation`).
    pub variant_classification: Option<String>,
    /// MAF `Consequence` — VEP Sequence Ontology terms. MAF packs multiple terms into
    /// one comma-separated cell, so this is a list.
    pub consequences: Vec<String>,
    /// MAF `Variant_Type` (`SNP`/`DEL`/`INS`/`DNP`/`TNP`/`ONP`).
    pub variant_type: Option<String>,
    pub hgvs_c: Option<String>,
    pub hgvs_p: Option<String>,
    /// `HGVSp_Short` — the `p.R1276*` form users actually search by.
    pub hgvs_p_short: Option<String>,
    pub protein_position: Option<String>,
    pub exon_number: Option<String>,
    pub impact: Option<String>,
    /// dbSNP ids (MAF `dbSNP_RS`; may be a comma/semicolon-separated list).
    pub dbsnp_rs: Vec<String>,
    pub gnomad_af: Option<f64>,
}

/// Read allele-depth columns; `vaf` is derived only when both numerator and denominator
/// are present and the denominator is non-zero.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Depth {
    pub t_depth: Option<u64>,
    pub t_ref_count: Option<u64>,
    pub t_alt_count: Option<u64>,
    pub n_depth: Option<u64>,
    pub n_ref_count: Option<u64>,
    pub n_alt_count: Option<u64>,
}

impl Depth {
    /// Tumor variant allele frequency, from `t_alt_count / t_depth`, falling back to
    /// `t_alt_count / (t_ref_count + t_alt_count)` when `t_depth` is absent.
    pub fn vaf(&self) -> Option<f64> {
        let alt = self.t_alt_count?;
        let total = match (self.t_depth, self.t_ref_count) {
            (Some(d), _) if d > 0 => d,
            (_, Some(r)) if r + alt > 0 => r + alt,
            _ => return None,
        };
        Some(alt as f64 / total as f64)
    }

    pub fn is_empty(&self) -> bool {
        *self == Depth::default()
    }
}

/// Split a MAF multi-value cell (VEP packs SO terms as `a,b`; dbSNP ids as `rs1;rs2`).
fn split_multi(v: Option<&str>) -> Vec<String> {
    v.map(|s| {
        s.split([',', ';', '&'])
            .map(str::trim)
            .filter(|t| !t.is_empty() && *t != "." && *t != "novel")
            .map(String::from)
            .collect()
    })
    .unwrap_or_default()
}

impl MafAnnotation {
    /// Read the annotation columns out of one MAF row.
    pub fn from_row(header: &MafHeader, row: &[&str]) -> Self {
        let s = |c: &str| header.get(row, c).map(String::from);
        Self {
            gene_symbol: s("Hugo_Symbol").filter(|g| g != "Unknown").or_else(|| s("SYMBOL")),
            entrez_gene_id: s("Entrez_Gene_Id").filter(|g| g != "0"),
            hgnc_id: s("HGNC_ID"),
            ensembl_gene_id: s("Gene"),
            transcript_id: s("Transcript_ID").or_else(|| s("Feature")),
            variant_classification: s("Variant_Classification"),
            consequences: split_multi(header.get(row, "Consequence")),
            variant_type: s("Variant_Type").or_else(|| s("VARIANT_CLASS")),
            hgvs_c: s("HGVSc"),
            hgvs_p: s("HGVSp"),
            hgvs_p_short: s("HGVSp_Short"),
            protein_position: s("Protein_position").or_else(|| s("Protein_Position")),
            exon_number: s("Exon_Number").or_else(|| s("EXON")),
            impact: s("IMPACT"),
            dbsnp_rs: split_multi(header.get(row, "dbSNP_RS")),
            gnomad_af: header
                .get(row, "gnomADg_AF")
                .or_else(|| header.get(row, "gnomAD_AF"))
                .and_then(|v| v.parse().ok()),
        }
    }
}

impl Depth {
    /// Read the depth columns out of one MAF row.
    pub fn from_row(header: &MafHeader, row: &[&str]) -> Self {
        let n = |c: &str| header.get(row, c).and_then(|v| v.parse().ok());
        Self {
            t_depth: n("t_depth"),
            t_ref_count: n("t_ref_count"),
            t_alt_count: n("t_alt_count"),
            n_depth: n("n_depth"),
            n_ref_count: n("n_ref_count"),
            n_alt_count: n("n_alt_count"),
        }
    }
}

/// One somatic observation of an allele in one tumor sample.
///
/// MAF has no `GT`, so there is no zygosity: a MAF row asserts *this sample carries this
/// allele*, with tumor allele depths as the quantitative evidence. The tumor/normal
/// barcode pair is what makes the call somatic, so both are carried.
#[derive(Debug, Clone)]
pub struct MafObservation {
    pub variant_id: String,
    pub tumor_sample: String,
    pub matched_normal: Option<String>,
    /// MAF `Chromosome`, as written in the file.
    pub contig: String,
    /// MAF 1-based `Start_Position` / `End_Position`.
    pub start: i64,
    pub end: Option<i64>,
    /// MAF-form alleles, `-` placeholders preserved, for round-tripping/debugging.
    pub reference_allele: String,
    pub alt_allele: String,
    pub assembly: Option<String>,
    pub reference_name: String,
    pub source: String,
    pub study_id: Option<String>,
    pub center: Option<String>,
    /// MAF `Mutation_Status`; NF-OSI MAFs leave it blank, so the CLI supplies a default
    /// (the whole cBioPortal `MUTATION_EXTENDED` profile is somatic).
    pub mutation_status: Option<String>,
    pub annotation: MafAnnotation,
    pub depth: Depth,
}

impl MafObservation {
    pub fn to_json(&self) -> Value {
        let mut m = Map::new();
        // NOTE: issue #95 models this class as `nf:VariantObservation`. The VCF path
        // still emits `VariantCall` (Beacon's vocabulary); unify before the RDF
        // projection lands.
        m.insert("type".into(), "VariantObservation".into());
        m.insert("variant".into(), self.variant_id.clone().into());
        m.insert("biosample".into(), self.tumor_sample.clone().into());
        m.insert("tumorSampleBarcode".into(), self.tumor_sample.clone().into());
        if let Some(n) = &self.matched_normal {
            m.insert("matchedNormalSampleBarcode".into(), n.clone().into());
        }
        m.insert("referenceName".into(), self.reference_name.clone().into());
        if let Some(a) = &self.assembly {
            m.insert("assemblyId".into(), a.clone().into());
        }
        m.insert("sourceContig".into(), self.contig.clone().into());
        m.insert("sourcePos".into(), self.start.into());
        if let Some(e) = self.end {
            m.insert("sourceEnd".into(), e.into());
        }
        m.insert("referenceBases".into(), self.reference_allele.clone().into());
        m.insert("alternateBases".into(), self.alt_allele.clone().into());
        m.insert("sourceFile".into(), self.source.clone().into());
        if let Some(s) = &self.study_id {
            m.insert("studyId".into(), s.clone().into());
        }
        if let Some(c) = &self.center {
            m.insert("center".into(), c.clone().into());
        }
        if let Some(s) = &self.mutation_status {
            m.insert("mutationStatus".into(), s.clone().into());
        }

        let a = &self.annotation;
        if let Some(g) = &a.ensembl_gene_id {
            m.insert("affectedGene".into(), g.clone().into());
        }
        if let Some(g) = &a.gene_symbol {
            m.insert("affectedGeneSymbol".into(), g.clone().into());
        }
        if let Some(g) = &a.entrez_gene_id {
            m.insert("entrezGeneId".into(), g.clone().into());
        }
        if let Some(g) = &a.hgnc_id {
            m.insert("hgncId".into(), g.clone().into());
        }
        if let Some(t) = &a.transcript_id {
            m.insert("transcriptId".into(), t.clone().into());
        }
        if let Some(v) = &a.variant_classification {
            m.insert("variantClassification".into(), v.clone().into());
        }
        if !a.consequences.is_empty() {
            m.insert("molecularConsequence".into(), a.consequences.clone().into());
        }
        if let Some(v) = &a.variant_type {
            m.insert("variantType".into(), v.clone().into());
        }
        if let Some(v) = &a.hgvs_c {
            m.insert("hgvsC".into(), v.clone().into());
        }
        if let Some(v) = &a.hgvs_p {
            m.insert("hgvsP".into(), v.clone().into());
        }
        // `aminoacidChange` mirrors the VCF path's key; HGVSp_Short is the searchable form.
        if let Some(v) = a.hgvs_p_short.as_ref().or(a.hgvs_p.as_ref()) {
            m.insert("aminoacidChange".into(), v.clone().into());
        }
        if let Some(v) = &a.protein_position {
            m.insert("proteinPosition".into(), v.clone().into());
        }
        if let Some(v) = &a.exon_number {
            m.insert("exonNumber".into(), v.clone().into());
        }
        if let Some(v) = &a.impact {
            m.insert("variantImpact".into(), v.clone().into());
        }
        if !a.dbsnp_rs.is_empty() {
            m.insert("dbsnpId".into(), a.dbsnp_rs.clone().into());
        }
        if let Some(v) = a.gnomad_af {
            m.insert("gnomadAlleleFrequency".into(), v.into());
        }

        let d = &self.depth;
        for (key, v) in [
            ("tumorDepth", d.t_depth),
            ("tumorRefCount", d.t_ref_count),
            ("tumorAltCount", d.t_alt_count),
            ("normalDepth", d.n_depth),
            ("normalRefCount", d.n_ref_count),
            ("normalAltCount", d.n_alt_count),
        ] {
            if let Some(v) = v {
                m.insert(key.into(), v.into());
            }
        }
        if let Some(vaf) = d.vaf() {
            m.insert("variantAlleleFrequency".into(), vaf.into());
        }
        Value::Object(m)
    }
}

/// A row that could not be given a VRS identity (unknown contig, assembly mismatch, or
/// non-ACGTN alleles). Per issue #95 these are **kept**, not dropped: a deterministic
/// local key from `{assembly}:{chrom}:{pos}:{ref}:{alt}` collapses duplicates across
/// studies, and the `unnormalized` flag lets the KG mark them `nf:unnormalizedVariant`.
#[derive(Debug, Clone)]
pub struct UnnormalizedVariant {
    pub key: String,
    pub assembly: String,
    pub contig: String,
    pub start: i64,
    pub reference_allele: String,
    pub alt_allele: String,
    pub reason: String,
}

impl UnnormalizedVariant {
    pub fn new(
        assembly: &str,
        contig: &str,
        start: i64,
        reference_allele: &str,
        alt_allele: &str,
        reason: impl Into<String>,
    ) -> Self {
        Self {
            key: format!("{assembly}:{contig}:{start}:{reference_allele}:{alt_allele}"),
            assembly: assembly.to_string(),
            contig: contig.to_string(),
            start,
            reference_allele: reference_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            reason: reason.into(),
        }
    }

    /// Local (non-VRS) id of the variant node. Observations of a rejected row reference
    /// this instead of a `ga4gh:VA.` id, so their sample/study provenance is not lost.
    pub fn id(&self) -> String {
        format!("nf:variant/{}", self.key)
    }

    pub fn to_json(&self) -> Value {
        let mut m = Map::new();
        m.insert("type".into(), "UnnormalizedVariant".into());
        m.insert("id".into(), self.id().into());
        m.insert("unnormalized".into(), true.into());
        m.insert("assemblyId".into(), self.assembly.clone().into());
        m.insert("sourceContig".into(), self.contig.clone().into());
        m.insert("sourcePos".into(), self.start.into());
        m.insert("referenceBases".into(), self.reference_allele.clone().into());
        m.insert("alternateBases".into(), self.alt_allele.clone().into());
        m.insert("reason".into(), self.reason.clone().into());
        Value::Object(m)
    }
}

/// Read a MAF, skipping `#` comment lines (cBioPortal MAFs carry
/// `#genome_nexus_version:` / `#isoform:` banners), and yield the header plus an
/// iterator-friendly line reader.
pub fn read_header<R: BufRead>(reader: &mut R) -> Result<MafHeader> {
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line).context("reading MAF header")?;
        if n == 0 {
            bail!("MAF has no header line");
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        return MafHeader::parse(&line);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "Hugo_Symbol\tEntrez_Gene_Id\tChromosome\tStart_Position\tEnd_Position\t\
        Variant_Type\tConsequence\tVariant_Classification\tReference_Allele\tTumor_Seq_Allele1\t\
        Tumor_Seq_Allele2\tTumor_Sample_Barcode\tNCBI_Build\tHGVSp_Short\tdbSNP_RS\tt_depth\t\
        t_ref_count\tt_alt_count";

    fn row(cells: &str) -> Vec<&str> {
        cells.split('\t').collect()
    }

    #[test]
    fn snv_interval_is_the_vcf_interval() {
        // A MAF SNP and the equivalent VCF record must land on the same interbase
        // interval — this is what makes the ids agree across formats.
        let (i, w) = maf_to_interbase(44908822, Some(44908822), "C", "T").unwrap();
        assert_eq!((i.start, i.end, i.alt.as_str()), (44908821, 44908822, "T"));
        assert!(w.is_none());
        assert!(!i.is_indel());
    }

    #[test]
    fn multi_nucleotide_substitution() {
        let (i, w) = maf_to_interbase(100, Some(102), "ACG", "TTT").unwrap();
        assert_eq!((i.start, i.end, i.alt.as_str()), (99, 102, "TTT"));
        assert!(w.is_none());
    }

    #[test]
    fn deletion_spans_the_deleted_bases() {
        // Real row: chr2 33539041-33539043 AAT -> "-"; VCF form anchors at 33539040.
        let (i, w) = maf_to_interbase(33539041, Some(33539043), "AAT", "-").unwrap();
        assert_eq!((i.start, i.end), (33539040, 33539043));
        assert_eq!(i.alt, "");
        assert_eq!(i.reference, "AAT");
        assert!(w.is_none());
        assert!(i.is_indel());
    }

    #[test]
    fn insertion_is_an_empty_interval_after_start() {
        // Real row: chr5 37325971-37325972 "-" -> A.
        let (i, w) = maf_to_interbase(37325971, Some(37325972), "-", "A").unwrap();
        assert_eq!((i.start, i.end), (37325971, 37325971));
        assert_eq!(i.alt, "A");
        assert_eq!(i.reference, "");
        assert!(w.is_none());
        assert!(i.is_indel());
    }

    #[test]
    fn disagreeing_end_position_is_reported_not_fatal() {
        let (i, w) = maf_to_interbase(100, Some(999), "ACG", "TTT").unwrap();
        // Reference_Allele wins.
        assert_eq!((i.start, i.end), (99, 102));
        assert!(w.unwrap().contains("expected 102"));
        // Absent End_Position is fine.
        assert!(maf_to_interbase(100, None, "A", "T").unwrap().1.is_none());
    }

    #[test]
    fn lowercase_alleles_are_uppercased() {
        let (i, _) = maf_to_interbase(100, Some(100), "c", "t").unwrap();
        assert_eq!((i.reference.as_str(), i.alt.as_str()), ("C", "T"));
    }

    #[test]
    fn effective_alt_prefers_allele2_then_allele1() {
        assert_eq!(effective_alt("C", "T", "C"), Some("T"));
        // Caller put the variant in allele1.
        assert_eq!(effective_alt("C", "C", "T"), Some("T"));
        // Germline-ref row: no variant to record.
        assert_eq!(effective_alt("C", "C", "C"), None);
        // Case-insensitive reference comparison.
        assert_eq!(effective_alt("C", "c", "T"), Some("T"));
    }

    #[test]
    fn assembly_aliases() {
        assert!(build_matches("GRCh38", "hg38"));
        assert!(build_matches("grch37", "hg19"));
        assert!(!build_matches("GRCh38", "GRCh37"));
        assert!(!build_matches("GRCh38", "hg19"));
    }

    #[test]
    fn header_requires_core_columns_and_accepts_variants() {
        // Older TCGA spelling is canonicalized.
        let h = MafHeader::parse(
            "Chromosome\tStart_position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode",
        )
        .unwrap();
        let r = row("17\t7577120\tC\tT\tTCGA-01");
        assert_eq!(h.get(&r, "Start_Position"), Some("7577120"));
        // Missing a required column is a hard error, not a silent empty run.
        let err = MafHeader::parse("Chromosome\tStart_Position").unwrap_err();
        assert!(err.to_string().contains("Reference_Allele"));
    }

    #[test]
    fn null_spellings_become_none() {
        let h = MafHeader::parse(HEADER).unwrap();
        let r = row(
            "NF1\t4763\t17\t31229110\t31229110\tSNP\tstop_gained\tNonsense_Mutation\tC\tC\tT\t\
             JH-2-001\tGRCh38\tp.R1276*\t.\t\t\t",
        );
        assert_eq!(h.get(&r, "HGVSp_Short"), Some("p.R1276*"));
        assert_eq!(h.get(&r, "dbSNP_RS"), None); // "."
        assert_eq!(h.get(&r, "t_depth"), None); // ""
    }

    #[test]
    fn annotation_reads_maf_columns() {
        let h = MafHeader::parse(HEADER).unwrap();
        let r = row(
            "NF1\t4763\t17\t31229110\t31229110\tSNP\t\
             splice_region_variant,intron_variant\tSplice_Region\tC\tC\tT\tJH-2-001\tGRCh38\t\
             p.R1276*\trs123;rs456\t60\t40\t20",
        );
        let a = MafAnnotation::from_row(&h, &r);
        assert_eq!(a.gene_symbol.as_deref(), Some("NF1"));
        assert_eq!(a.entrez_gene_id.as_deref(), Some("4763"));
        assert_eq!(a.variant_classification.as_deref(), Some("Splice_Region"));
        // MAF packs several SO terms in one cell; each must survive for SSSOM mapping.
        assert_eq!(a.consequences, vec!["splice_region_variant", "intron_variant"]);
        assert_eq!(a.hgvs_p_short.as_deref(), Some("p.R1276*"));
        assert_eq!(a.dbsnp_rs, vec!["rs123", "rs456"]);

        let d = Depth::from_row(&h, &r);
        assert_eq!(d.t_alt_count, Some(20));
        assert_eq!(d.vaf(), Some(20.0 / 60.0));
    }

    #[test]
    fn unknown_gene_symbol_is_dropped() {
        let h = MafHeader::parse(HEADER).unwrap();
        let r = row(
            "Unknown\t0\t17\t31229110\t31229110\tSNP\t\tIGR\tC\tC\tT\tJH-2-001\tGRCh38\t\t\t\t\t",
        );
        let a = MafAnnotation::from_row(&h, &r);
        assert_eq!(a.gene_symbol, None);
        assert_eq!(a.entrez_gene_id, None); // Entrez 0 is MAF's "no gene"
    }

    #[test]
    fn vaf_falls_back_to_ref_plus_alt() {
        let d = Depth {
            t_ref_count: Some(30),
            t_alt_count: Some(10),
            ..Default::default()
        };
        assert_eq!(d.vaf(), Some(0.25));
        // No alt count → no VAF, and never a divide-by-zero.
        assert_eq!(
            Depth {
                t_depth: Some(0),
                t_alt_count: Some(5),
                ..Default::default()
            }
            .vaf(),
            None
        );
        assert_eq!(Depth::default().vaf(), None);
    }

    #[test]
    fn header_skips_cbioportal_comment_banners() {
        let mut r = std::io::Cursor::new(
            "#genome_nexus_version: 1.0.2\n#isoform: mskcc\n\
             Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n",
        );
        let h = read_header(&mut r).unwrap();
        assert!(h.has("Tumor_Sample_Barcode"));
    }

    #[test]
    fn unnormalized_key_is_deterministic() {
        let a = UnnormalizedVariant::new("GRCh38", "1", 100, "A", "-", "unknown contig");
        let b = UnnormalizedVariant::new("GRCh38", "1", 100, "A", "-", "different reason");
        assert_eq!(a.key, b.key);
        assert_eq!(a.key, "GRCh38:1:100:A:-");
    }
}
