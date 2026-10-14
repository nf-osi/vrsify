//! vrsify — assign GA4GH VRS 2.0 identity to variants from cohort files.
//!
//! Reads VCF ([`vcf`]) or MAF ([`maf`]) and emits GA4GH VRS 2.0 alleles with computed
//! identifiers plus per-sample observation records, for the NF Beacon KG pipeline. The
//! identity engine ([`digest`], [`vrs`], [`normalize`], [`refget`]) is format-agnostic;
//! only the front ends are format-specific, and **the same variant gets the same
//! `ga4gh:VA.` id through either one**.

pub mod digest;
pub mod maf;
pub mod normalize;
pub mod refget;
pub mod unnormalized;
pub mod vcf;
pub mod vrs;
