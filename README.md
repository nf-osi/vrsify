# vrsify

**Give every variant in a VCF or MAF a stable, global GA4GH VRS identifier.**

`vrsify` reads a VCF or a cBioPortal-style MAF (`data_mutations.txt`) and writes two
NDJSON streams:

1. **alleles** — one [GA4GH VRS 2.0](https://vrs.ga4gh.org) object per distinct variant,
   with its computed `ga4gh:VA.…` id;
2. **observations** — one record per sample (or per tumor sample, for MAF) pointing at
   that id, carrying the context-dependent details: genotype/zygosity, read depths, VAF,
   study and file provenance, and any VEP/snpEff annotation.

Because a VRS id is essentially a digest of the variant itself (assembly sequence, position, 
alternate bases), the *same* edit always hashes to the *same* id, no matter which file
or which format it came from. That makes it a great join key with no additional coordinate 
or representation bookkeeping needed.

```
$ vrsify --vcf cohort.vcf.gz --seqmap seqmap.tsv --reference GRCh38.fa \
         --out-alleles alleles.ndjson --out-observations obs.ndjson

# alleles.ndjson  (context-free: exactly id + location + state + type)
{"id":"ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
 "location":{"start":44908821,"end":44908822,
             "sequenceReference":{"refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                                  "type":"SequenceReference"},
             "id":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","type":"SequenceLocation"},
 "state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}

# obs.ndjson  (context-full: everything that is about a sample, not the variant)
{"variant":"ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt","biosample":"S1",
 "zygosity":"heterozygous","referenceName":"19","assemblyId":"GRCh38",
 "sourcePos":44908822,"sourceContig":"chr19","referenceBases":"C","alternateBases":"T",
 "sourceFile":"cohort.vcf.gz","type":"VariantObservation"}
```

## Why

The reference implementation [**ga4gh/vrs-python**](https://github.com/ga4gh/vrs-python) 
is used for VRS definition and verifying correctness; **`vrsify` is a Rust implementation
that additionally covers cBioPortal maf format, is faster, and has several other advantages**. 

`vrsify` was written for these use case:

- **No Python environment, no SeqRepo, no UTA.** vrs-python resolves sequences through a
  data proxy backed by a local SeqRepo instance (tens of GB) or a REST service. `vrsify`
  is a single static binary that needs only a reference FASTA and a small TSV mapping
  each contig to its refget accession — it will *generate* that TSV from the FASTA
  for you. Nothing to install, nothing to keep running, works offline and in CI.
- **MAF is a first-class input.** Much public somatic data (cBioPortal studies, TCGA) is
  distributed as MAF, not VCF. vrs-python's annotator targets VCF, so the usual route is
  MAF → VCF → VRS, which needs a reference lookup just to re-create the anchor base that
  MAF dropped. `vrsify` projects MAF straight into VRS coordinates, and a variant
  ingested from a MAF gets *the same* `ga4gh:VA.` id as the same variant from a VCF. That
  equivalence is asserted in the tests.
- **Cohort-scale throughput.** Streaming, allele-deduplicating, single pass: a real
  23,741-row cBioPortal MAF converts in about half a second.
- **Output shaped for loading, not annotating.** Rather than writing VRS ids back into
  INFO fields, `vrsify` emits the two-layer split a graph or data warehouse actually wants:
  the context-free allele (shared across every cohort that sees it) and the per-sample
  observation. Sample, genotype, study, and annotation data are *never* placed on the
  allele.

While this was built for the NF knowledge-graph pipeline, where cBioPortal MAF
studies and VCF-derived calls have to join on variant identity, this should be generally reusable.

**Use vrs-python instead if** you need HGVS/SPDI/gnomAD-string translation, transcript
projection, or the parts of the VRS model `vrsify` doesn't cover here (see
[Scope](#scope-and-limitations)).

## Install

Requires a Rust toolchain ([rustup](https://rustup.rs)).

```sh
cargo build --release      # binary at ./target/release/vrsify
cargo test                 # includes the GA4GH golden-fixture gate
```

## Quick start

### 1. Make a seqmap from your reference

VRS identifies a sequence by its **refget accession** (`SQ.` + a digest of the sequence
itself), not by a contig name, so each contig in your input has to be mapped to one:

```sh
vrsify seqmap --fasta GRCh38.fa --assembly GRCh38 --out seqmap.tsv
```

`seqmap.tsv` is a plain TSV — `contig`, `refgetAccession`, and optionally `assemblyId`
and `referenceName` — so you can also hand-write it or pull accessions from elsewhere:

```
chr19	SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl	GRCh38	19
```

`assemblyId` and `referenceName` are copied onto observations as the familiar
human-readable coordinates; `refgetAccession` is what identity is computed from. A
`chr` prefix mismatch between your input and your seqmap (`19` vs `chr19`) is resolved
automatically.

### 2. Convert a VCF

```sh
vrsify \
  --vcf input.vcf.gz \            # plain or bgzipped, detected by content
  --seqmap seqmap.tsv \
  --reference GRCh38.fa \         # optional, but see "Indels" below
  --out-alleles alleles.ndjson \
  --out-observations obs.ndjson \
  --source "syn12345/input.vcf"   # free-text provenance, recorded on each observation
```

Split multi-allelic sites and left-align upstream (`bcftools norm -m- -f ref.fa`) as
usual; `vrsify` then applies the VRS-specific normalization on top. Alleles are
case-insensitive per the VCF spec, so `a>t` and `A>T` get the same id.

`zygosity` is read from `GT`: `homozygous`, `heterozygous`, `hemizygous` (a haploid
call such as `GT=1`), or `unknown` when part of the genotype is uncalled (`GT=1/.` is
diploid with one unknown allele, so the copy count is not knowable).

### 3. Convert a MAF

```sh
vrsify maf \
  --maf data_mutations.txt \      # `#` banner lines are skipped
  --seqmap seqmap.tsv \
  --reference GRCh38.fa \         # required for correct indel ids
  --out-alleles alleles.ndjson \
  --out-observations obs.ndjson \
  --study-id nst_nfosi_ntap \
  --source "cbioportal:nst_nfosi_ntap/data_mutations.txt"
```

MAF observations carry the tumor/normal barcode pair, allele depths and derived VAF, the
study id, and MAF's own annotation columns (`Variant_Classification`, `Consequence` as a
**list** of SO terms, `HGVSc`/`HGVSp`/`HGVSp_Short`, `Transcript_ID`, `Gene`, `HGNC_ID`,
`dbSNP_RS`, gnomAD AF). MAF has no `GT`, so there is no zygosity.

## Indels, and why `--reference` matters

Substitutions are byte-exact without a reference: the common REF/ALT prefix and suffix
are trimmed either way (that part of VRS normalization needs no sequence lookup), so even
a padded spelling like `GT>GA` reaches the canonical id of the underlying `T>A`. Indels
are not: the same insertion or
deletion inside a homopolymer or tandem repeat can be written at several positions, and
VRS resolves that by **fully justifying** the variant — rolling it as far left and as far
right as the reference allows — before hashing. That requires the reference sequence.

With `--reference`, equivalent representations collapse to one id, matching
vrs-python and ClinVar. Without it, indel alleles are emitted with
`"fullyJustified": false` (a flag in the output; not part of the digest) and the count is
printed as a warning. Substitutions are unaffected either way. Only about 4% of a typical
MAF's rows are indels, but those ids will not join anything unless you pass the FASTA.

The FASTA is currently loaded per-contig into memory; a whole-genome run needs roughly
3 GB of RAM.

When `--reference` is supplied, both converters require each converted contig to be
present in the FASTA, verify its sequence digest against the seqmap accession, and
check that the variant interval is within bounds and its REF bases agree with the
FASTA (case-insensitively). A mismatch stops conversion with an error; supplying a
reference never silently falls back to an unnormalized projection. Without
`--reference`, both VCF and MAF indels carry the `fullyJustified: false` marker and a
warning.

## How MAF rows become VRS coordinates

MAF already stores indels trimmed, using `-` placeholders, so the projection into VRS's
interbase (0-based, half-open) coordinates is a pure coordinate shift — no anchor-base
lookup, the step a MAF → VCF conversion needs a reference for:

| MAF row                                           | VRS interbase | state   |
|---------------------------------------------------|---------------|---------|
| SNP `Start=100 End=102 ref=ACG alt=TTT`           | `[99, 102)`   | `"TTT"` |
| DEL `Start=100 End=102 ref=ACG alt=-`             | `[99, 102)`   | `""`    |
| INS `Start=100 End=101 ref=- alt=TT`              | `[100, 100)`  | `"TT"`  |

`Reference_Allele` (not `End_Position`) is authoritative for the span, since it is what
the alleles are built from; a disagreeing `End_Position` is counted and sampled into the
run summary rather than failing the row.

## Data-quality behaviour

`vrsify` is deliberately loud and drops nothing you didn't ask it to drop.

- **Missing alleles are not deletions.** Missing `Reference_Allele` or
  `Tumor_Seq_Allele2` values (blank, `.`, or `NA`) stop MAF conversion with a row-specific
  error. Only an explicit `-` means an empty allele. Reference validation failures
  also stop conversion. Output files from a failed run may contain earlier rows and
  must not be treated as complete.
- **Assembly is part of identity.** A row whose `NCBI_Build` disagrees with the seqmap's
  `assemblyId` never gets a VRS id (handled assemblies are `hg38`/`GRCh38`,
  `hg19`/`GRCh37` and `hg18`/`NCBI36`, each with any `.p13`-style GRC patch suffix — a
  patch release adds scaffolds without moving primary-assembly coordinates, so
  `GRCh37.p13` and `GRCh37` are the same coordinate system and must not read as a
  mismatch). Real studies do mix builds — one we ingested had 29 GRCh37 rows inside an
  otherwise GRCh38 study, which would otherwise have been hashed against GRCh38
  accessions and given plausible, wrong ids.
- **Nothing is discarded.** Records that can't be given an identity — unknown contig,
  assembly mismatch (MAF), non-nucleotide alleles — are still written to the alleles
  stream as `{"type":"UnnormalizedVariant","unnormalized":true, …,"reason":…}` with a
  deterministic `{assembly}:{chrom}:{pos}:{ref}:{alt}` key, and counted in the run
  summary. Their observations are written too, referencing that `nf:variant/…` id, so
  the sample, study and source of a rejected record survive even though its VRS
  identity does not — two samples carrying the same unknown-contig variant still yield
  two observations. The node and observations use the source/CLI assembly when present,
  otherwise the resolved seqmap assembly, with `unknown` used only when neither can
  identify it. **This holds for both front ends**, and `--strict` turns it into a hard
  failure in either. On the VCF path, pass `--assembly` so unknown-contig keys are
  something better than `unknown:…` (a seqmap `assemblyId` always wins over it).
- **No silent filtering.** `--min-tumor-alt-count` is **off** by default. It exists
  because real MAFs contain rows with `t_alt_count = 0` — no read in the tumor supports
  the allele (28% of rows in one study we ingested). Recording those as "this specimen
  carries this variant" would be wrong, so pass `--min-tumor-alt-count 1` when that
  applies; nothing is filtered unless you ask.
- **Structural variants are out of scope, explicitly.** Symbolic and structural VCF
  ALTs (`<DEL>`, breakends) and the `*` spanning-deletion placeholder are the one thing
  that is skipped rather than kept: a `<DEL>` is defined by `INFO/END` or `SVLEN`, which
  a `{ref}:{alt}`-keyed node would silently drop, and `*` is not an allele of its own.
  They are counted in the run summary. VCF records with no `POS` are likewise counted
  and skipped — there is no locus to key a variant on.

## Correctness

VRS ids are worthless unless they are byte-identical to everyone else's. 
`cargo test` checks:

- the `sha512t24u` digest primitive against GA4GH's `functions.yaml` vectors;
- every `Allele`, `SequenceLocation`, `CopyNumberCount`, `CopyNumberChange`, and
  `Adjacency` case in GA4GH's vendored `vrs@2.0` `models.yaml` golden suite —
  serialization, digest, and computed id, **byte-exact**;
- **cross-format identity**: a MAF SNP row reproduces the GA4GH golden id, and a deletion
  and an insertion written in MAF's trimmed form and in VCF's anchored form collapse to
  one id;
- fully-justified normalization (left- and right-shifted representations of one indel
  must converge) and `seqmap` generation, on synthetic references;
- missing MAF alleles, reference/seqmap identity and REF mismatches, interval bounds,
  identity-allele normalization, and warnings/markers for VCF indels without a reference;
- reference-free trimming (padded substitutions reach the canonical id with or without
  `--reference`, in both formats), and that identity alleles (`REF == ALT`) reach one id
  either way;
- that both front ends keep unknown-contig and non-ACGTN records as
  `UnnormalizedVariant` nodes with their observations, honour `--strict`, and label
  observations with the same `type`;
- GRC patch-release builds (`GRCh37.p13`) matching their base assembly, and per-ALT
  `CSQ` selection + percent-decoding on multiallelic records.

Two checks need multi-GB reference files and so are not part of `cargo test`, but were
run against real data: `vrsify seqmap` over GRCh38 chr19 reproduces GA4GH's canonical
accession for that contig (`SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl`), and 13 real chr19
variants (SNV, unique indel, single- and multi-unit STR indels) produce ids identical to
`ga4gh.vrs` 2.3.3 — 13/13.

Reproduce end-to-end from the bundled examples:

```sh
vrsify --vcf examples/sample.vcf --seqmap examples/seqmap.tsv \
  --out-alleles examples/alleles.ndjson --out-observations examples/obs.ndjson \
  --source "syn12345/test.vcf"
grep -q 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt' examples/alleles.ndjson \
  && echo "rs7412 VA id matches the GA4GH golden fixture"

vrsify maf --maf examples/sample.maf --seqmap examples/seqmap.tsv \
  --out-alleles examples/maf_alleles.ndjson --out-observations examples/maf_obs.ndjson \
  --study-id nst_nfosi_ntap --source "cbioportal:nst_nfosi_ntap/data_mutations.txt"
grep -q 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt' examples/maf_alleles.ndjson \
  && echo "the MAF path reproduces the same id"
```

## Scope and limitations

Supported VRS types: `Allele` with `SequenceLocation` / `SequenceReference` and
`LiteralSequenceExpression` / `ReferenceLengthExpression` / `LengthExpression` states,
plus `CopyNumberCount`, `CopyNumberChange`, and `Adjacency` in the library (the
CNV/SV types are id-gated against the golden suite; the VCF and MAF front ends currently
emit alleles only). `CisPhasedBlock`, `DerivativeMolecule`, and `Terminus` are not
modelled.

Not implemented: live refget/SeqRepo service lookups (everything is derived from your
FASTA), HGVS or SPDI translation, transcript-level projection.

VEP `CSQ` and snpEff `ANN` INFO annotations, when present, are carried onto each
observation as `affectedGene` / `affectedGeneSymbol` / `aminoacidChange` /
`molecularConsequence` (never onto the context-free allele). A `CSQ`/`ANN` value holds
one entry per (allele, transcript), so entries are first filtered to the ALT being
observed — by `ALLELE_NUM` when VEP wrote it, otherwise by the allele column, allowing
for VEP's minimal (`-`-padded) indel spelling — and the first survivor is taken, which
is the most-severe one. That matters on multiallelic records and on records split by
`bcftools norm -m-`, which leaves the full original `CSQ` on every split line. Field
values are percent-decoded per VCF 4.3 (`p.Cys130%3D` → `p.Cys130=`). Annotation is
advisory metadata; it never affects the VRS id.

## License

MIT — see [LICENSE](LICENSE).
