//! Reference-based "fully justified" VRS normalization.
//!
//! Ports the algorithm used by `vrs-python` (`ga4gh.vrs.normalize`, `_normalize` for the
//! `LiteralSequenceExpression` case): trim the common prefix/suffix of REF vs ALT, then,
//! for pure insertions/deletions, *bidirectionally justify* the indel across any repeat
//! by rolling left and right as far as the flanking reference allows. This yields a
//! canonical interbase interval + replacement sequence so that, e.g., all equivalent
//! representations of an indel in a homopolymer/tandem repeat collapse to one VRS id —
//! matching what `bcftools norm` + `vrs-python` produce.
//!
//! Coarse left-alignment and multiallelic splitting are expected upstream
//! (`bcftools norm -m- -f REF`); this is the VRS-specific final step.
//!
//! ## Validation status
//! Unit tests here exercise the algorithm over **small synthetic references** (a
//! homopolymer and a tandem-repeat region) that demonstrate bidirectional justification.
//! Separately, the emitted ids were cross-checked against `ga4gh.vrs` 2.3.3 over 13 real
//! GRCh38 chr19 variants (SNV, unique indel, single- and multi-unit STR indels): all 13
//! `ga4gh:VA.` ids matched byte-for-byte. That check needs the multi-GB reference FASTA,
//! so it is not part of `cargo test`.

use anyhow::{ensure, Result};

/// Read-only access to a reference sequence, in interbase (0-based, half-open)
/// coordinates. `len` is the total sequence length; `get(a, b)` returns the residues in
/// `[a, b)` (uppercased ASCII). Callers validate intervals before requesting residues.
pub trait Reference {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Return residues in interbase `[start, end)`. Callers guarantee `start <= end` and
    /// `end <= len()`.
    fn get(&self, start: usize, end: usize) -> Vec<u8>;

    /// Sequence identity used to verify the seqmap before minting an allele id.
    /// Large reference providers should cache this value.
    fn refget_accession(&self) -> String {
        crate::refget::refget_accession(&self.get(0, self.len()))
    }
}

/// An in-memory reference sequence, handy for tests and small contigs.
pub struct InMemoryReference {
    seq: Vec<u8>,
    refget: String,
}

impl InMemoryReference {
    pub fn new(seq: impl Into<Vec<u8>>) -> Self {
        let mut seq = seq.into();
        seq.make_ascii_uppercase();
        let refget = format!("SQ.{}", crate::digest::sha512t24u(&seq));
        Self { seq, refget }
    }
}

impl Reference for InMemoryReference {
    fn len(&self) -> usize {
        self.seq.len()
    }
    fn get(&self, start: usize, end: usize) -> Vec<u8> {
        self.seq[start.min(self.seq.len())..end.min(self.seq.len())].to_vec()
    }

    fn refget_accession(&self) -> String {
        self.refget.clone()
    }
}

/// Check that the reference belongs to the seqmap and agrees with the input edit.
pub fn validate_reference(
    reference: &dyn Reference,
    start: usize,
    end: usize,
    expected_ref: &[u8],
    expected_accession: &str,
) -> Result<()> {
    validate_interval(reference, start, end)?;
    let accession = reference.refget_accession();
    ensure!(
        accession == expected_accession,
        "reference accession {accession} does not match seqmap accession {expected_accession}"
    );
    ensure!(
        reference.get(start, end).eq_ignore_ascii_case(expected_ref),
        "REF mismatch at interbase [{start}, {end}): input REF '{}' disagrees with reference FASTA",
        String::from_utf8_lossy(expected_ref)
    );
    Ok(())
}

fn validate_interval(reference: &dyn Reference, start: usize, end: usize) -> Result<()> {
    ensure!(
        start <= end && end <= reference.len(),
        "interval [{start}, {end}) is outside reference bounds (length {})",
        reference.len()
    );
    Ok(())
}

/// A set of reference contigs loaded from a FASTA, keyed by contig name. Used at ingest
/// time to justify indels against the reference.
///
/// NOTE: this loads whole contigs into memory. That is fine for the small synthetic
/// fixtures used in tests and for modest references, but a full human genome (~3 GB)
/// should instead use `noodles-fasta`'s indexed (`.fai`) windowed `query` around each
/// variant. Left as a TODO since large references are unavailable in this environment.
pub struct FastaReferences {
    contigs: std::collections::HashMap<String, InMemoryReference>,
}

impl FastaReferences {
    /// Load all contigs from a (plain or bgzipped) FASTA file.
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Self> {
        use noodles_fasta as fasta;
        let mut reader = fasta::io::reader::Builder.build_from_path(path)?;
        let mut contigs = std::collections::HashMap::new();
        for result in reader.records() {
            let record = result?;
            let name = String::from_utf8_lossy(record.name()).to_string();
            contigs.insert(
                name,
                InMemoryReference::new(record.sequence().as_ref().to_vec()),
            );
        }
        Ok(Self { contigs })
    }

    pub fn contig(&self, name: &str) -> Option<&InMemoryReference> {
        self.contigs.get(name)
    }
}

/// The canonical VRS state chosen for a normalized allele, matching vrs-python's
/// `_normalize_allele` state selection (see module docs and `normalize`).
///
/// * `Literal` — a `LiteralSequenceExpression{sequence: alt}` (substitutions/SNVs and
///   pure insertions that do not expand a reference repeat).
/// * `ReferenceLengthExpression` — an indel over a non-empty reference region, carrying
///   the fully-justified `length` (= `alt.len()`) and `repeat_subunit_length`. The `alt`
///   bytes are still available on [`NormalizedAllele::alt`] for the RLE `sequence` field
///   (which is excluded from the VRS digest by `vrs.rs`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizedState {
    Literal,
    ReferenceLengthExpression {
        length: usize,
        repeat_subunit_length: usize,
    },
}

/// A normalized allele: an interbase interval on the reference, the replacement (alt)
/// sequence for that interval, and the canonical VRS state kind chosen for it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedAllele {
    pub start: usize,
    pub end: usize,
    pub alt: Vec<u8>,
    pub state: NormalizedState,
}

/// Fully-justified VRS normalization of a REF→ALT edit at interbase `[start, end)`.
///
/// Callers with source REF bases must first call [`validate_reference`]. Invalid
/// intervals and empty no-op edits return errors. Identity alleles retain their
/// original span as a ReferenceLengthExpression.
pub fn normalize(
    reference: &dyn Reference,
    start: usize,
    end: usize,
    alt: &[u8],
) -> Result<NormalizedAllele> {
    validate_interval(reference, start, end)?;
    let ref_bases = reference.get(start, end);
    let alt = alt.to_ascii_uppercase();

    // Identity alleles have no moving indel unit; do not enter repeat expansion.
    if ref_bases == alt {
        ensure!(!alt.is_empty(), "cannot normalize an empty no-op edit");
        return Ok(NormalizedAllele {
            start,
            end,
            state: NormalizedState::ReferenceLengthExpression {
                length: alt.len(),
                repeat_subunit_length: alt.len(),
            },
            alt,
        });
    }

    // 1. Trim common suffix, then common prefix. Track how the interval shrinks.
    let (pfx, sfx) = trim_common(&ref_bases, &alt);
    let ref_seq = ref_bases[pfx..ref_bases.len() - sfx].to_vec();
    let alt_seq = alt[pfx..alt.len() - sfx].to_vec();
    let new_start = start + pfx;
    let new_end = end - sfx;

    // 2. If both sides are non-empty after trimming, it's a substitution / MNV — VRS
    //    does not roll these; return the trimmed form as a LiteralSequenceExpression
    //    (vrs-python step 2.b).
    if !ref_seq.is_empty() && !alt_seq.is_empty() {
        return Ok(NormalizedAllele {
            start: new_start,
            end: new_end,
            alt: alt_seq,
            state: NormalizedState::Literal,
        });
    }

    // `seed_length` is the length of the trimmed indel *before* EXPAND rolling — the
    // deleted reference length for a deletion, the inserted length for an insertion.
    // vrs-python uses it directly as the RLE `repeatSubunitLength` for deletions, and as
    // the pool of candidate factors for the insertion cycle search. (`seed_length =
    // len_trimmed_ref or len_trimmed_alt`.)
    let seed_length = if ref_seq.is_empty() { alt_seq.len() } else { ref_seq.len() };

    // 3. Pure insertion (ref empty) or deletion (alt empty): bidirectionally justify by
    //    EXPANDing the interval over the full ambiguous region (vrs-python / bioutils
    //    EXPAND mode). The "unit" is the moving sequence — the inserted bases for an
    //    insertion, the deleted reference bases for a deletion.
    Ok(justify_expand(
        reference, new_start, new_end, &ref_seq, &alt_seq, seed_length,
    ))
}

/// The reference-free part of VRS normalization: the lengths of the common suffix and
/// common prefix (trimmed in that order, matching vrs-python) shared by a REF/ALT pair.
/// Returns `(prefix_len, suffix_len)`; the trims never overlap. The no-reference
/// projections use this too, so a padded substitution (`GT>GA`) reaches the same trimmed
/// interval — and therefore the same VRS id — with or without a reference FASTA.
pub fn trim_common(reference: &[u8], alt: &[u8]) -> (usize, usize) {
    let mut sfx = 0;
    while sfx < reference.len().min(alt.len())
        && reference[reference.len() - 1 - sfx] == alt[alt.len() - 1 - sfx]
    {
        sfx += 1;
    }
    let (r, a) = (&reference[..reference.len() - sfx], &alt[..alt.len() - sfx]);
    let mut pfx = 0;
    while pfx < r.len().min(a.len()) && r[pfx] == a[pfx] {
        pfx += 1;
    }
    (pfx, sfx)
}

/// bioutils/vrs-python EXPAND-mode justification. `del` is the (trimmed) deleted
/// reference segment, `ins` the (trimmed) inserted segment; exactly one is empty for a
/// pure indel. Rolls the moving unit left and right as far as the flanking reference
/// permits, then expands the interval to cover the whole ambiguous region and prepends /
/// appends the uncovered reference bases to the allele.
fn justify_expand(
    reference: &dyn Reference,
    start: usize,
    end: usize,
    del: &[u8],
    ins: &[u8],
    seed_length: usize,
) -> NormalizedAllele {
    // The repeat unit that can circularly permute is the non-empty side.
    let unit = if del.is_empty() { ins } else { del };
    let ulen = unit.len();
    debug_assert!(ulen > 0);

    // roll_left: max distance d such that, extending left, the cyclically-indexed unit
    // matches the reference. Compare unit[-(k+1) % ulen] to ref[start-1-k].
    let mut ldist = 0usize;
    while start - ldist > 0 {
        let ref_base = reference.get(start - ldist - 1, start - ldist);
        let unit_base = unit[(ulen - 1) - (ldist % ulen)];
        if ref_base.first().copied() == Some(unit_base) {
            ldist += 1;
        } else {
            break;
        }
    }
    // roll_right: compare unit[k % ulen] to ref[end + k].
    let mut rdist = 0usize;
    while end + rdist < reference.len() {
        let ref_base = reference.get(end + rdist, end + rdist + 1);
        let unit_base = unit[rdist % ulen];
        if ref_base.first().copied() == Some(unit_base) {
            rdist += 1;
        } else {
            break;
        }
    }

    let new_start = start - ldist;
    let new_end = end + rdist;
    let lseq = reference.get(new_start, start);
    let rseq = reference.get(end, new_end);

    // Build the fully-justified alt (`extended_alt`) and the reference region it spans
    // (`extended_ref = reference[new_start, new_end)`).
    let alt = if del.is_empty() {
        // Insertion: the reference interval is the fully-justified ambiguous flank
        // [new_start, new_end); the alt is that flank with the inserted unit added.
        let mut a = lseq;
        a.extend_from_slice(ins);
        a.extend_from_slice(&rseq);
        a
    } else {
        // Deletion: the reference interval covers the ambiguous region including the
        // deleted unit; the alt is the flanking sequence with the unit removed.
        let mut a = lseq;
        a.extend_from_slice(&rseq);
        a
    };
    let extended_ref = reference.get(new_start, new_end);
    let state = choose_state(&extended_ref, &alt, seed_length);
    NormalizedAllele {
        start: new_start,
        end: new_end,
        alt,
        state,
    }
}

/// Select the canonical VRS state for a fully-justified pure indel, mirroring
/// vrs-python's `_normalize_allele` steps 5.a–5.c exactly.
///
/// * 5.a — empty reference region (a pure insertion that did not expand a repeat) →
///   `Literal`.
/// * 5.b — deletion (`alt` shorter than the reference region) →
///   `ReferenceLengthExpression{length: alt.len(), repeatSubunitLength: seed_length}`.
///   (The subunit length is the *trimmed deleted length*, NOT the minimal repeat period:
///   deleting two copies of a k-mer yields `repeatSubunitLength = 2k`.)
/// * 5.c — insertion (`alt` longer than the reference region) → the largest factor `f` of
///   `seed_length` with `f <= extended_ref.len()` for which the inserted tail is a valid
///   cyclic continuation of `extended_ref[extended_ref.len()-f..]` gives
///   `ReferenceLengthExpression{length: alt.len(), repeatSubunitLength: f}`; if no factor
///   qualifies, `Literal`.
fn choose_state(extended_ref: &[u8], alt: &[u8], seed_length: usize) -> NormalizedState {
    // 5.a
    if extended_ref.is_empty() {
        return NormalizedState::Literal;
    }
    // 5.b: deletion.
    if alt.len() < extended_ref.len() {
        return NormalizedState::ReferenceLengthExpression {
            length: alt.len(),
            repeat_subunit_length: seed_length,
        };
    }
    // 5.c: insertion. Try factors of seed_length in descending order.
    if alt.len() > extended_ref.len() {
        for cycle_length in factors_desc(seed_length) {
            if cycle_length > extended_ref.len() {
                continue;
            }
            let cycle_start = extended_ref.len() - cycle_length;
            if is_valid_cycle(cycle_start, extended_ref, alt) {
                return NormalizedState::ReferenceLengthExpression {
                    length: alt.len(),
                    repeat_subunit_length: cycle_length,
                };
            }
        }
        // 5.c.3: no valid cycle → literal.
        return NormalizedState::Literal;
    }
    // Equal lengths cannot occur for a pure indel post-expansion; treat as literal.
    NormalizedState::Literal
}

/// All factors of `n` in descending order (matches vrs-python `_factor_gen`, which yields
/// largest-first). Returns an empty vector for `n == 0`.
fn factors_desc(n: usize) -> Vec<usize> {
    let mut upper = Vec::new();
    let mut lower = Vec::new();
    let mut i = 1;
    while i * i <= n {
        if n.is_multiple_of(i) {
            upper.push(n / i);
            if n / i != i {
                lower.push(i);
            }
        }
        i += 1;
    }
    lower.reverse();
    upper.extend(lower);
    upper
}

/// vrs-python `_is_valid_cycle`: the portion of `target` past `template.len()` must be a
/// cyclic continuation of `template[template_start..]`.
fn is_valid_cycle(template_start: usize, template: &[u8], target: &[u8]) -> bool {
    let cycle = &template[template_start..];
    if cycle.is_empty() {
        return template.len() >= target.len();
    }
    for (k, &ch) in target[template.len()..].iter().enumerate() {
        if ch != cycle[k % cycle.len()] {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_allele_keeps_original_span_without_repeat_expansion() {
        let r = InMemoryReference::new("CAAAAG");
        let n = normalize(&r, 2, 4, b"AA").unwrap();
        assert_eq!((n.start, n.end), (2, 4));
        assert_eq!(n.alt, b"AA");
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression {
                length: 2,
                repeat_subunit_length: 2,
            }
        );
    }

    #[test]
    fn invalid_intervals_and_empty_no_ops_return_errors() {
        let r = InMemoryReference::new("ACGT");
        for (start, end, alt) in [
            (3, 2, b"A".as_slice()),
            (3, 5, b"A"),
            (5, 5, b"A"),
            (2, 2, b""),
        ] {
            assert!(normalize(&r, start, end, alt).is_err());
        }
    }

    #[test]
    fn snv_unchanged() {
        // ref ...A[C]G... , C>T at interbase [3,4)
        let r = InMemoryReference::new("AAACGTT");
        let n = normalize(&r, 3, 4, b"T").unwrap();
        assert_eq!(
            n,
            NormalizedAllele {
                start: 3,
                end: 4,
                alt: b"T".to_vec(),
                state: NormalizedState::Literal
            }
        );
    }

    #[test]
    fn trims_common_bases_then_expands() {
        // VCF-style padded deletion: REF=CTT ALT=CT at [0,3) → delete one T.
        // reference: C T T A. Trimming yields a 1-base T deletion, then EXPAND covers
        // the whole ambiguous TT run [1,3); alt keeps one T.
        let r = InMemoryReference::new("CTTA");
        let n = normalize(&r, 0, 3, b"CT").unwrap();
        assert_eq!((n.start, n.end), (1, 3));
        assert_eq!(n.alt, b"T".to_vec()); // one of the two T's remains
        assert_eq!(r.get(n.start, n.end), b"TT");
        // Deletion over a non-empty ref region → RLE, length=1, repeatSubunitLength=seed(1).
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression { length: 1, repeat_subunit_length: 1 }
        );
    }

    #[test]
    fn deletion_fully_justified_in_homopolymer() {
        // reference: G AAAAA C  (5 A's at [1,6)). Delete one A at [3,4).
        let r = InMemoryReference::new("GAAAAAC");
        let n = normalize(&r, 3, 4, b"").unwrap();
        // EXPAND covers the entire A-run; alt holds 4 A's (one deleted).
        assert_eq!((n.start, n.end), (1, 6));
        assert_eq!(n.alt, b"AAAA".to_vec());
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression { length: 4, repeat_subunit_length: 1 }
        );
    }

    #[test]
    fn insertion_fully_justified_in_tandem_repeat() {
        // reference: C ATAT G. Insert one "AT" unit at [1,1) into the AT tandem repeat.
        let r = InMemoryReference::new("CATATG");
        let n = normalize(&r, 1, 1, b"AT").unwrap();
        // EXPAND covers the AT repeat [1,5); alt is three AT units (one inserted).
        assert_eq!((n.start, n.end), (1, 5));
        assert_eq!(n.alt, b"ATATAT".to_vec());
        // Insertion expanding a tandem repeat → RLE; seed=2 (one "AT" unit), period 2.
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression { length: 6, repeat_subunit_length: 2 }
        );
    }

    #[test]
    fn deletion_of_two_repeat_units_uses_seed_length_as_subunit() {
        // reference: G (AT)x5 C. Delete two "AT" units → seed_length = 4 (the trimmed
        // deleted length). vrs-python emits repeatSubunitLength = 4, NOT the period 2.
        let r = InMemoryReference::new("GATATATATATC");
        // Delete "ATAT" at [1,5) (two units).
        let n = normalize(&r, 1, 5, b"").unwrap();
        assert_eq!((n.start, n.end), (1, 11)); // whole AT run [1,11)
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression { length: 6, repeat_subunit_length: 4 }
        );
    }

    #[test]
    fn insertion_of_two_repeat_units_uses_full_seed_as_subunit() {
        // reference: G (AT)x5 C. Insert two "AT" units inside the run → seed_length = 4,
        // largest valid cycle factor of 4 is 4 → repeatSubunitLength = 4.
        let r = InMemoryReference::new("GATATATATATC");
        let n = normalize(&r, 3, 3, b"ATAT").unwrap();
        assert_eq!((n.start, n.end), (1, 11));
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression { length: 14, repeat_subunit_length: 4 }
        );
    }

    #[test]
    fn trim_common_trims_suffix_then_prefix() {
        assert_eq!(trim_common(b"GT", b"GA"), (1, 0));
        assert_eq!(trim_common(b"TG", b"AG"), (0, 1));
        assert_eq!(trim_common(b"ACGT", b"AGGT"), (1, 2));
        assert_eq!(trim_common(b"A", b"T"), (0, 0));
        // Padded indels: the suffix is consumed first, matching `normalize`.
        assert_eq!(trim_common(b"AA", b"A"), (0, 1));
        assert_eq!(trim_common(b"A", b"AA"), (0, 1));
        // Identity is consumed entirely by the suffix pass; the trims never overlap.
        assert_eq!(trim_common(b"AT", b"AT"), (0, 2));
        assert_eq!(trim_common(b"", b""), (0, 0));
    }

    #[test]
    fn factors_desc_matches_reference() {
        assert_eq!(factors_desc(1), vec![1]);
        assert_eq!(factors_desc(2), vec![2, 1]);
        assert_eq!(factors_desc(4), vec![4, 2, 1]);
        assert_eq!(factors_desc(6), vec![6, 3, 2, 1]);
        assert_eq!(factors_desc(12), vec![12, 6, 4, 3, 2, 1]);
    }

    #[test]
    fn insertion_rolls_into_flanking_base() {
        let r = InMemoryReference::new("ACGTACGT");
        // Insert "TTT" at [4,4): rolls left one base into the T at index 3.
        let n = normalize(&r, 4, 4, b"TTT").unwrap();
        assert_eq!((n.start, n.end), (3, 4));
        assert_eq!(n.alt, b"TTTT".to_vec());
        // Insertion into a 1-base flank: extended_ref="T" (non-empty) so RLE; the largest
        // valid cycle factor of seed 3 that fits the 1-base ref is 1.
        assert_eq!(
            n.state,
            NormalizedState::ReferenceLengthExpression { length: 4, repeat_subunit_length: 1 }
        );
    }

    #[test]
    fn insertion_no_repeat_unchanged() {
        let r = InMemoryReference::new("ACGTACGT");
        // Insert "CCC" at [1,1): C at index 1 differs; no left roll (ref[0]='A'),
        // no right roll (ref[1]='C' vs unit[0]='C' → actually matches!). Use a unit that
        // cannot roll at all: insert "GGG" at [4,4) — ref[3]='T', ref[4]='A', no roll.
        let n = normalize(&r, 4, 4, b"GGG").unwrap();
        assert_eq!((n.start, n.end), (4, 4));
        assert_eq!(n.alt, b"GGG".to_vec());
        // Pure insertion with empty reference region → LiteralSequenceExpression.
        assert_eq!(n.state, NormalizedState::Literal);
    }
}
