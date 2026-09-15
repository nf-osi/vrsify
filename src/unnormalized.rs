//! Variants that could not be given a VRS identity, kept rather than dropped.
//!
//! Both front ends ([`crate::vcf`] and [`crate::maf`]) reject some records from VRS
//! identity: the contig is not in the seqmap, the source assembly disagrees with the
//! seqmap's, or the alleles are not plain nucleotides. Dropping those rows would
//! silently lose the sample/study provenance attached to them, so instead they become
//! an `UnnormalizedVariant` node with a deterministic local id, and their observations
//! are emitted pointing at that id.
//!
//! ## The id is local, so its namespace is the caller's to choose
//!
//! A normalized allele gets a `ga4gh:VA.` id that means the same thing to everyone. An
//! unnormalized one cannot: all that is left is a deterministic key over the source
//! coordinates, which is only unique *within* whoever minted it. `vrsify` therefore
//! refuses to guess a prefix for it — the caller supplies one (`--variant-id-prefix`),
//! and a run that would need one without having been given one fails rather than
//! stamping someone else's data with a namespace it has no claim to.

use serde_json::{Map, Value};

/// A row that could not be given a VRS identity (unknown contig, assembly mismatch, or
/// non-ACGTN alleles). Per issue #95 these are **kept**, not dropped: a deterministic
/// local key from `{assembly}:{chrom}:{pos}:{ref}:{alt}` collapses duplicates across
/// studies, and the `unnormalized` flag lets the KG mark them `nf:unnormalizedVariant`.
#[derive(Debug, Clone)]
pub struct UnnormalizedVariant {
    /// Caller-supplied id namespace, used verbatim (see [`UnnormalizedVariant::id`]).
    pub id_prefix: String,
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
        id_prefix: &str,
        assembly: &str,
        contig: &str,
        start: i64,
        reference_allele: &str,
        alt_allele: &str,
        reason: impl Into<String>,
    ) -> Self {
        Self {
            id_prefix: id_prefix.to_string(),
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
    ///
    /// The prefix is concatenated verbatim, so the caller controls the separator:
    /// `nf:variant/` yields `nf:variant/GRCh38:1:100:A:-`.
    pub fn id(&self) -> String {
        format!("{}{}", self.id_prefix, self.key)
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
        let a = UnnormalizedVariant::new("nf:variant/", "GRCh38", "1", 100, "A", "-", "unknown contig");
        let b = UnnormalizedVariant::new("nf:variant/", "GRCh38", "1", 100, "A", "-", "different reason");
        // The key is over the source coordinates only: the reason a row was rejected is
        // metadata, so two studies rejecting the same row collapse to one node.
        assert_eq!(a.key, b.key);
        assert_eq!(a.key, "GRCh38:1:100:A:-");
        assert_eq!(a.id(), "nf:variant/GRCh38:1:100:A:-");
    }

    #[test]
    fn the_id_prefix_is_used_verbatim() {
        let v = |prefix| {
            UnnormalizedVariant::new(prefix, "GRCh38", "1", 100, "A", "T", "unknown contig").id()
        };
        // Whatever separator the caller wants — including none at all.
        assert_eq!(v("ex:var/"), "ex:var/GRCh38:1:100:A:T");
        assert_eq!(v("https://example.org/variant/"), "https://example.org/variant/GRCh38:1:100:A:T");
        assert_eq!(v(""), "GRCh38:1:100:A:T");
    }
}
