//! Variants that could not be given a VRS identity, kept rather than dropped.
//!
//! Both front ends ([`crate::vcf`] and [`crate::maf`]) reject some records from VRS
//! identity: the contig is not in the seqmap, the source assembly disagrees with the
//! seqmap's, or the alleles are not plain nucleotides. Dropping those rows would
//! silently lose the sample/study provenance attached to them, so instead they become
//! an `UnnormalizedVariant` node with a deterministic local id, and their observations
//! are emitted pointing at that id.

use serde_json::{Map, Value};

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unnormalized_key_is_deterministic() {
        let a = UnnormalizedVariant::new("GRCh38", "1", 100, "A", "-", "unknown contig");
        let b = UnnormalizedVariant::new("GRCh38", "1", 100, "A", "-", "different reason");
        assert_eq!(a.key, b.key);
        assert_eq!(a.key, "GRCh38:1:100:A:-");
        assert_eq!(a.id(), "nf:variant/GRCh38:1:100:A:-");
    }
}
