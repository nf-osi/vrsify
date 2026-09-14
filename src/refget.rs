//! Refget accession computation + seqmap generation.
//!
//! The GA4GH refget / VRS `SequenceReference.refgetAccession` for a sequence is
//! `SQ.` + `sha512t24u(uppercased_sequence_bytes)` — the same `sha512t24u` primitive
//! used for VRS object digests, applied to the residue string (uppercased, no
//! whitespace). This lets us derive a `--seqmap` directly from a reference FASTA instead
//! of hand-supplying it, and to detect the assembly per-file.
//!
//! ## Validation status
//! The digest *algorithm* (uppercase → SHA-512 → first 24 bytes → base64url) is the
//! GA4GH-specified refget identity and is unit-tested here on a small synthetic
//! sequence, and the `sha512t24u` primitive is gated byte-exact against GA4GH's
//! `functions.yaml`. Separately, running `vrsify seqmap` over the real GRCh38 chr19
//! FASTA reproduced GA4GH's canonical accession for that contig
//! (`SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl`); that check needs the multi-GB FASTA, so it
//! is not part of `cargo test`.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_fasta as fasta;

use crate::digest::sha512t24u;

/// Compute the refget accession (`SQ.<digest>`) for a residue sequence. The sequence is
/// uppercased before hashing, per refget/VRS convention.
pub fn refget_accession(sequence: &[u8]) -> String {
    let mut upper = sequence.to_vec();
    upper.make_ascii_uppercase();
    format!("SQ.{}", sha512t24u(&upper))
}

/// One contig's derived seqmap row.
#[derive(Debug, Clone)]
pub struct SeqmapRow {
    pub contig: String,
    pub refget: String,
    pub assembly: Option<String>,
    pub reference_name: String,
    pub length: usize,
}

/// Normalize a contig name to a Beacon `referenceName` (strip a leading `chr`).
pub fn reference_name_for(contig: &str) -> String {
    contig.strip_prefix("chr").unwrap_or(contig).to_string()
}

/// Read a reference FASTA and compute a seqmap row per contig. `assembly` is applied to
/// every row (per-file assembly is a caller concern; supply `None` to leave it blank).
///
/// Streams sequence bytes so it does not require the whole FASTA in memory at once (per
/// record it holds one contig's sequence — unavoidable for the digest).
pub fn seqmap_from_fasta<P: AsRef<Path>>(
    path: P,
    assembly: Option<&str>,
) -> Result<Vec<SeqmapRow>> {
    let path = path.as_ref();
    let mut reader = fasta::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("opening reference FASTA {path:?}"))?;

    let mut rows = Vec::new();
    for result in reader.records() {
        let record = result.context("reading FASTA record")?;
        let contig = String::from_utf8_lossy(record.name()).to_string();
        let seq = record.sequence().as_ref();
        rows.push(SeqmapRow {
            refget: refget_accession(seq),
            length: seq.len(),
            reference_name: reference_name_for(&contig),
            assembly: assembly.map(String::from),
            contig,
        });
    }
    Ok(rows)
}

/// Write seqmap rows in the TSV format `load_seqmap` reads:
/// `contig <tab> refgetAccession <tab> assemblyId <tab> referenceName`.
pub fn write_seqmap<W: Write>(mut out: W, rows: &[SeqmapRow]) -> Result<()> {
    writeln!(out, "# contig\trefgetAccession\tassemblyId\treferenceName")?;
    for r in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            r.contig,
            r.refget,
            r.assembly.as_deref().unwrap_or(""),
            r.reference_name,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refget_matches_primitive() {
        // sha512t24u("ACGT") is a known functions.yaml vector; the SQ accession is the
        // primitive prefixed with "SQ.".
        assert_eq!(refget_accession(b"ACGT"), "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
        // Case-insensitive: lowercase input hashes to the same accession.
        assert_eq!(refget_accession(b"acgt"), refget_accession(b"ACGT"));
    }

    #[test]
    fn reference_name_strips_chr() {
        assert_eq!(reference_name_for("chr19"), "19");
        assert_eq!(reference_name_for("19"), "19");
        assert_eq!(reference_name_for("chrX"), "X");
    }
}
