# vrsify

Generate stable GA4GH VRS identifiers for variants in VCF and MAF files.

`vrsify` reads a VCF or a cBioPortal-style MAF (`data_mutations.txt`) and writes two
NDJSON streams:

1. **alleles** — one [GA4GH VRS 2.0](https://vrs.ga4gh.org) object per distinct variant,
   with its computed `ga4gh:VA.…` identifier. This layer is the standard.
2. **observations** — one `VariantObservation` per sample, containing genotype or tumor
   measurements, provenance, and VEP/snpEff annotations. This is a `vrsify`-specific
   format because VRS does not model sample-level observations.

Because a VRS id is essentially a digest of the variant itself (assembly sequence, position, 
alternate bases), the *same* edit always hashes to the *same* id, no matter which file
or which format it came from. That makes it a great join key with no additional coordinate 
or representation bookkeeping needed.

```
$ vrsify --vcf cohort.vcf.gz --seqmap seqmap.tsv --reference GRCh38.fa \
         --out-alleles alleles.ndjson --out-observations obs.ndjson \
         --strict

# alleles.ndjson (context-free VRS object with exactly id + location + state + type)
{"id":"ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
 "location":{"start":44908821,"end":44908822,
             "sequenceReference":{"refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                                  "type":"SequenceReference"},
             "id":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","type":"SequenceLocation"},
 "state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}

# obs.ndjson (sample-specific observation)
{"variant":"ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt","biosample":"S1",
 "zygosity":"heterozygous","referenceName":"19","assemblyId":"GRCh38",
 "sourcePos":44908822,"sourceContig":"chr19","referenceBases":"C","alternateBases":"T",
 "sourceFile":"cohort.vcf.gz","type":"VariantObservation"}
```

## Design goals

[ga4gh/vrs-python](https://github.com/ga4gh/vrs-python) is the reference implementation
used to validate `vrsify` output. `vrsify` is a Rust converter for these goals and use cases:

- **Minimal runtime dependencies.** vrs-python resolves sequences through a
  data proxy backed by a local SeqRepo instance (tens of GB) or a REST service. `vrsify`
  is a single static binary that needs only a reference FASTA and a small TSV mapping
  each contig to its refget accession — it will *generate* that TSV from the FASTA. 
  Nothing to install, nothing to keep running, works offline and in CI.
- **MAF is a first-class input.** Much public somatic data (cBioPortal studies, TCGA) is
  distributed as MAF, not VCF. vrs-python's annotator targets VCF, so the usual route is
  MAF → VCF → VRS, which needs a reference lookup just to re-create the anchor base that
  MAF dropped. `vrsify` projects MAF straight into VRS coordinates, and a variant 
  ingested from a MAF gets *the same* `ga4gh:VA.` id as the same variant from a VCF.
- **Cohort-scale throughput.** Streaming conversion where input is processed in a single pass, 
  with duplicate alleles removed from the allele stream. A typical cBioPortal MAF 
  converts in about half a second.
- **Loading-oriented output.** Context-free alleles and sample-specific observations are
  emitted separately for optimized loading into graphs or data warehouses.

Use vrs-python for HGVS, SPDI, or gnomAD-string translation; transcript projection; or
VRS types outside the [supported scope](#scope-and-limitations).

## Install

Requires a Rust toolchain ([rustup](https://rustup.rs)).

```sh
cargo build --release  # binary: ./target/release/vrsify
cargo test             # includes GA4GH conformance fixtures
```

## Quick start

### 1. Make a seqmap from a reference

VRS identifies a sequence by its **refget accession** (`SQ.` + a digest of the sequence
itself), not by a contig name. Each input contig must therefore be mapped to an
accession:

```sh
vrsify seqmap --fasta GRCh38.fa --assembly GRCh38 --out seqmap.tsv
```

`seqmap.tsv` is a plain TSV — `contig`, `refgetAccession`, and optionally `assemblyId`
and `referenceName`. It can also be created manually with accessions from another
source:

```
chr19	SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl	GRCh38	19
```

`assemblyId` and `referenceName` are copied to observations as human-readable
coordinates; VRS identities use `refgetAccession`. A `chr` prefix mismatch between the
input and seqmap (`19` vs `chr19`) is resolved automatically.

### 2. Convert a VCF

```sh
# Input may be plain text or bgzipped; compression is detected from the content.
vrsify \
  --vcf input.vcf.gz \
  --seqmap seqmap.tsv \
  --reference GRCh38.fa \
  --out-alleles alleles.ndjson \
  --out-observations obs.ndjson \
  --source "example-study/input.vcf" \
  --strict
```

`--reference` is optional for substitutions but required for fully justified indel
identifiers; see [Indels and reference sequences](#indels-and-reference-sequences).
`--source` is free-text provenance recorded on each observation. `--strict` is required
unless `--variant-id-prefix` is provided; see
[Data-quality behavior](#data-quality-behavior).

Split multiallelic sites and left-align them before conversion, for example with
`bcftools norm -m- -f ref.fa`. `vrsify` then applies VRS normalization. VCF alleles are
case-insensitive, so `a>t` and `A>T` produce the same identifier.

`zygosity` is read from `GT`: `homozygous`, `heterozygous`, `hemizygous` (a haploid
call such as `GT=1`), or `unknown` when part of the genotype is uncalled (`GT=1/.` is
diploid with one unknown allele, so the copy count cannot be determined).

### 3. Convert a MAF

```sh
vrsify maf \
  --maf data_mutations.txt \ # `#` banner lines are skipped
  --seqmap seqmap.tsv \
  --reference GRCh38.fa \
  --out-alleles alleles.ndjson \
  --out-observations obs.ndjson \
  --study-id example-study \
  --source "cbioportal:example-study/data_mutations.txt" \
  --strict # or provide --variant-id-prefix
```

MAF observations include the tumor and normal barcodes, allele depths, derived VAF,
study identifier, and annotation columns (`Variant_Classification`, `Consequence` as a
list of SO terms, `HGVSc`/`HGVSp`/`HGVSp_Short`, `Transcript_ID`, `Gene`, `HGNC_ID`,
`dbSNP_RS`, and gnomAD AF). MAF does not define `GT`, so zygosity is not emitted.

## The two output layers

`vrsify` separates context-free variants from sample-specific observations.

**`alleles.ndjson`** contains one GA4GH VRS `Allele` per distinct variant, identified by
a stable `ga4gh:VA.` digest. These identifiers can be used to join equivalent variants
across files, formats, and other VRS-compatible sources. The stream may also contain:

- `"fullyJustified": false` — an extra key on an otherwise conformant `Allele`, marking
  an indel normalized without `--reference`. The key is excluded from the digest and is
  omitted when a reference is provided.
- `"type":"UnnormalizedVariant"` — a record that could not be given a VRS identity at
  all. Its source information is retained; with `--strict`, the conversion fails
  instead.

With `--reference --strict`, this stream contains only VRS `Allele` objects.

**`obs.ndjson`** contains one `VariantObservation` per sample and allele. Because VRS
does not model sample-level calls, this stream is a `vrsify`-specific loading format. It
includes zygosity or tumor depth and VAF, provenance, and VEP/snpEff annotations. The
`variant` field references an identifier in `alleles.ndjson`, including for
unnormalized variants. 

## Indels and reference sequences

Substitution normalization trims common REF/ALT prefixes and suffixes without a sequence
lookup. For example, `GT>GA` produces the canonical identifier for the underlying
`T>A` substitution.

Indels in repeated sequence can have several equivalent representations. VRS fully
justifies an indel by extending it left and right against the reference before computing
the identifier. This requires `--reference`.

With a reference, equivalent indel representations produce the same identifier as
vrs-python. Without a reference, indel alleles include `"fullyJustified": false`, which
is excluded from the digest, and the run reports their count as a warning. Substitution
identifiers are unaffected.

The FASTA is currently loaded per-contig into memory; a whole-genome run needs roughly
3 GB of RAM.

When `--reference` is supplied, both converters require each converted contig to be
present in the FASTA, verify its sequence digest against the seqmap accession, and
check that the variant interval is within bounds and its REF bases agree with the
FASTA (case-insensitively). A mismatch stops conversion with an error; supplying a
reference does not fall back to an unnormalized projection.

## How MAF rows become VRS coordinates

MAF stores trimmed indels with `-` placeholders. Conversion to VRS interbase (0-based,
half-open) coordinates requires a coordinate shift but no anchor-base lookup:

| MAF row                                           | VRS interbase | state   |
|---------------------------------------------------|---------------|---------|
| SUB `Start=100 End=102 ref=ACG alt=TTT`           | `[99, 102)`   | `"TTT"` |
| DEL `Start=100 End=102 ref=ACG alt=-`             | `[99, 102)`   | `""`    |
| INS `Start=100 End=101 ref=- alt=TT`              | `[100, 100)`  | `"TT"`  |

`Reference_Allele` determines the span. If `End_Position` disagrees, the row is still
processed and the discrepancy is included in the run summary.

## Data-quality behavior

- **Missing MAF alleles cause an error.** Blank, `.`, or `NA` values in
  `Reference_Allele` or `Tumor_Seq_Allele2` stop conversion with a row-specific error.
  Only `-` represents an empty allele. Reference validation errors also stop conversion.
  Output from a failed run may be incomplete.
- **Assembly is part of variant identity.** A MAF `NCBI_Build` value that conflicts with
  the seqmap `assemblyId` cannot produce a VRS identifier. Supported aliases are
  `hg38`/`GRCh38`, `hg19`/`GRCh37`, and `hg18`/`NCBI36`. GRC patch suffixes such as
  `.p13` are accepted because they do not change primary-assembly coordinates.
- **Unnormalized records can be retained.** Unknown contigs, MAF assembly mismatches,
  and non-nucleotide alleles are written as `UnnormalizedVariant` records when
  `--variant-id-prefix` is provided. Each record receives a deterministic
  `{assembly}:{chrom}:{pos}:{ref}:{alt}` local key, is counted in the run summary, and
  retains its observations and provenance. `--strict` instead stops conversion. For
  VCF input, `--assembly` supplies the assembly for unknown-contig keys; a seqmap
  `assemblyId` takes precedence.
- **Local identifiers require a namespace.** Each run must provide either
  `--variant-id-prefix` or `--strict`. The prefix is used verbatim, including its
  separator. For example, `--variant-id-prefix example:variant/` produces
  `example:variant/GRCh38:chr1:100:A:T`. This requirement is validated before the
  input is opened.
- **Tumor-alt filtering is optional.** `--min-tumor-alt-count` is disabled by default.
  Set `--min-tumor-alt-count 1` to exclude MAF rows with no supporting tumor reads.
- **Unsupported VCF records are skipped.** Symbolic and structural ALTs (`<DEL>` and
  breakends), the `*` spanning-deletion placeholder, and records without `POS` are
  counted in the run summary but not emitted. These records cannot be represented by
  the local key used for unnormalized small variants.

## Correctness

`cargo test` validates:

- the `sha512t24u` digest primitive against GA4GH's `functions.yaml` vectors;
- every `Allele`, `SequenceLocation`, `CopyNumberCount`, `CopyNumberChange`, and
  `Adjacency` case in the GA4GH VRS 2.0 `models.yaml` conformance fixtures, including
  serialization, digest, and computed identifier;
- cross-format identity: a MAF SNP row reproduces the GA4GH reference identifier, and
  deletion and insertion records written in MAF's trimmed form and VCF's anchored form
  collapse to one identifier;
- fully-justified normalization (left- and right-shifted representations of one indel
  must converge) and `seqmap` generation, on synthetic references;
- missing MAF alleles, reference/seqmap identity and REF mismatches, interval bounds,
  identity-allele normalization, and warnings/markers for VCF indels without a reference;
- reference-free trimming (padded substitutions reach the canonical identifier with or
  without `--reference`, in both formats), and that identity alleles (`REF == ALT`)
  reach one identifier either way;
- both converters retain unknown-contig and non-ACGTN records as
  `UnnormalizedVariant` nodes with their observations, honor `--strict`, use the same
  observation `type`, and require an explicit namespace for local identifiers before
  reading the input;
- GRC patch-release builds (`GRCh37.p13`) matching their base assembly, and per-ALT
  `CSQ` selection + percent-decoding on multiallelic records.

Two validations require multi-GB reference files and are not part of `cargo test`.
`vrsify seqmap` over GRCh38 chromosome 19 reproduces the canonical GA4GH accession
`SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl`. Thirteen chromosome 19 variants covering an SNV,
a unique indel, and single- and multi-unit STR indels produced identifiers identical to
`ga4gh.vrs` 2.3.3.

The bundled examples support an end-to-end check:

```sh
vrsify --vcf examples/sample.vcf --seqmap examples/seqmap.tsv \
  --out-alleles examples/alleles.ndjson --out-observations examples/obs.ndjson \
  --source "example-study/test.vcf" --strict
grep -q 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt' examples/alleles.ndjson \
  && echo "rs7412 identifier matches the GA4GH fixture"

vrsify maf --maf examples/sample.maf --seqmap examples/seqmap.tsv \
  --out-alleles examples/maf_alleles.ndjson --out-observations examples/maf_obs.ndjson \
  --study-id example-study --source "cbioportal:example-study/data_mutations.txt" \
  --strict
grep -q 'ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt' examples/maf_alleles.ndjson \
  && echo "MAF identifier matches the VCF identifier"
```

## Scope and limitations

Supported VRS types: `Allele` with `SequenceLocation` / `SequenceReference` and
`LiteralSequenceExpression` / `ReferenceLengthExpression` / `LengthExpression` states,
plus `CopyNumberCount`, `CopyNumberChange`, and `Adjacency` in the library (the
CNV/SV identifiers are validated against the same conformance fixtures; the VCF and MAF
converters emit alleles only). `CisPhasedBlock`, `DerivativeMolecule`, and `Terminus`
are not modeled.

Not implemented: live refget/SeqRepo service lookups, HGVS or SPDI translation, and
transcript-level projection. Sequence data is read from the local FASTA.

VEP `CSQ` and snpEff `ANN` INFO annotations, when present, are copied to each
observation as `affectedGene` / `affectedGeneSymbol` / `aminoacidChange` /
`molecularConsequence`, not to the context-free allele. Because a `CSQ`/`ANN` value has
one entry per allele and transcript, entries are filtered to the observed ALT using
`ALLELE_NUM` when available or the allele column otherwise. VEP's minimal,
`-`-padded indel representation is supported. The first remaining entry, representing
the most severe consequence, is used. This selection applies to multiallelic records
and records split by `bcftools norm -m-`, which retain the original `CSQ` entries. Field
values are percent-decoded according to VCF 4.3 (`p.Cys130%3D` → `p.Cys130=`).
Annotations do not affect the VRS identifier.

## License

MIT — see [LICENSE](LICENSE).
