//! vrsify CLI: VCF/MAF → VRS 2.0 alleles (NDJSON) + observations (NDJSON).
//!
//! I/O uses `noodles-vcf`: transparently reads plain-text and bgzipped (`.vcf.gz`)
//! VCF, with spec-conformant header/sample parsing. The `seqmap` subcommand derives
//! `SQ.` refget accessions from a reference FASTA. The `maf` subcommand ingests
//! cBioPortal-style MAF (`data_mutations.txt`) through the same VRS engine, emitting
//! somatic observations instead of genotype calls.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use noodles_vcf as vcf;

use vrsify::maf::{
    allele_bases, build_maf_allele, build_matches, effective_alt, maf_to_interbase, read_header,
    Depth, MafAnnotation, MafObservation,
};
use vrsify::normalize::FastaReferences;
use vrsify::refget::{seqmap_from_fasta, write_seqmap};
use vrsify::unnormalized::UnnormalizedVariant;
use vrsify::vcf::{
    build_allele, build_allele_normalized, csq_format_from_description, extract_annotation,
    is_concrete_alt, load_seqmap, resolve_contig, zygosity_for, Observation, SeqMap,
};

#[derive(Parser, Debug)]
#[command(name = "vrsify", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    convert: ConvertArgs,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Convert a VCF to VRS alleles + observations (default when no subcommand given).
    Convert(ConvertArgs),
    /// Convert a MAF (cBioPortal `data_mutations.txt`) to VRS alleles + somatic observations.
    Maf(MafArgs),
    /// Derive a `--seqmap` TSV of `SQ.` refget accessions from a reference FASTA.
    Seqmap(SeqmapArgs),
}

#[derive(Parser, Debug, Default)]
struct ConvertArgs {
    /// Input VCF, plain or bgzipped (`.vcf.gz`). Detected by content.
    #[arg(long)]
    vcf: Option<PathBuf>,

    /// TSV seqmap: contig<tab>refgetAccession[<tab>assemblyId[<tab>referenceName]].
    #[arg(long)]
    seqmap: Option<PathBuf>,

    /// Output NDJSON of unique VRS alleles.
    #[arg(long)]
    out_alleles: Option<PathBuf>,

    /// Output NDJSON of per-sample observations.
    #[arg(long)]
    out_observations: Option<PathBuf>,

    /// Reference FASTA (optional). When supplied, indels are fully-justified normalized
    /// against the validated reference; otherwise indels are marked fullyJustified=false.
    #[arg(long)]
    reference: Option<PathBuf>,

    /// Identifier/URI for the source VCF, recorded on each observation.
    #[arg(long, default_value = "")]
    source: String,

    /// Assembly to record when the seqmap does not declare one — in particular on
    /// unnormalized records for contigs the seqmap does not contain, whose
    /// `{assembly}:{chrom}:{pos}:{ref}:{alt}` key would otherwise say `unknown`.
    /// A seqmap `assemblyId` always wins over this.
    #[arg(long)]
    assembly: Option<String>,

    /// Fail instead of emitting unnormalized records (unknown contig / non-ACGTN
    /// alleles).
    #[arg(long)]
    strict: bool,

    /// Id namespace for records that cannot be given a VRS id, used verbatim (include
    /// the separator: `nf:variant/` yields `nf:variant/GRCh38:chr1:100:A:T`). Such an
    /// id is local — it means something only inside the namespace that minted it, so
    /// there is no correct default and one must be named. Required unless `--strict`,
    /// which rejects those records outright and so needs no local id.
    #[arg(long)]
    variant_id_prefix: Option<String>,
}

#[derive(Parser, Debug)]
struct SeqmapArgs {
    /// Reference FASTA (plain or bgzipped).
    #[arg(long)]
    fasta: PathBuf,

    /// Assembly id to stamp on every row (e.g. GRCh38). Optional.
    #[arg(long)]
    assembly: Option<String>,

    /// Output TSV path; defaults to stdout.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Parser, Debug)]
struct MafArgs {
    /// Input MAF (tab-separated; `#` banner lines are skipped).
    #[arg(long)]
    maf: PathBuf,

    /// TSV seqmap: contig<tab>refgetAccession[<tab>assemblyId[<tab>referenceName]].
    /// `chr`-prefix differences between the MAF and the seqmap are resolved automatically.
    #[arg(long)]
    seqmap: PathBuf,

    /// Output NDJSON of unique variant nodes (VRS alleles, plus any unnormalized rows).
    #[arg(long)]
    out_alleles: PathBuf,

    /// Output NDJSON of per-tumor-sample somatic observations.
    #[arg(long)]
    out_observations: PathBuf,

    /// Reference FASTA. REQUIRED for indel ids to be fully justified (and therefore to
    /// match vrs-python / ClinVar / the VCF path). Without it, indels are emitted with
    /// `"fullyJustified": false` and counted.
    #[arg(long)]
    reference: Option<PathBuf>,

    /// Identifier/URI for the source MAF, recorded on each observation.
    #[arg(long, default_value = "")]
    source: String,

    /// cBioPortal study id (e.g. `nst_nfosi_ntap`), recorded on each observation.
    #[arg(long)]
    study_id: Option<String>,

    /// Assembly to assume when the MAF has no `NCBI_Build` column. When the MAF *does*
    /// carry a build and the seqmap declares one, a mismatch sends the row to the
    /// unnormalized stream rather than minting a wrong-assembly VRS id.
    #[arg(long)]
    assembly: Option<String>,

    /// Value for `mutationStatus` when the MAF's own `Mutation_Status` is blank.
    /// cBioPortal `MUTATION_EXTENDED` profiles are somatic by construction.
    #[arg(long, default_value = "Somatic")]
    mutation_status: String,

    /// Drop rows whose `t_alt_count` is below this. Off by default so nothing is
    /// filtered silently. 28% of `nfib_ctf_biobank_2025` rows have `t_alt_count = 0`
    /// (no read in the tumor supports the allele) — recording those as "this specimen
    /// carries this variant" would be wrong, so set `--min-tumor-alt-count 1` (or
    /// higher) for that study. Rows with no `t_alt_count` column at all are kept.
    #[arg(long, default_value_t = 0)]
    min_tumor_alt_count: u64,

    /// Fail instead of emitting unnormalized rows (unknown contig / assembly mismatch /
    /// non-ACGTN alleles).
    #[arg(long)]
    strict: bool,

    /// Id namespace for rows that cannot be given a VRS id, used verbatim (include the
    /// separator: `nf:variant/` yields `nf:variant/GRCh38:chr1:100:A:T`). Such an id is
    /// local — it means something only inside the namespace that minted it, so there is
    /// no correct default and one must be named. Required unless `--strict`, which
    /// rejects those rows outright and so needs no local id.
    #[arg(long)]
    variant_id_prefix: Option<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Some(Command::Seqmap(args)) => run_seqmap(args),
        Some(Command::Maf(args)) => run_maf(args),
        Some(Command::Convert(args)) => run_convert(args),
        None => run_convert(cli.convert),
    }
}

fn run_seqmap(args: SeqmapArgs) -> Result<()> {
    let rows = seqmap_from_fasta(&args.fasta, args.assembly.as_deref())?;
    match &args.out {
        Some(p) => {
            let f = File::create(p).with_context(|| format!("creating {p:?}"))?;
            write_seqmap(BufWriter::new(f), &rows)?;
        }
        None => write_seqmap(std::io::stdout().lock(), &rows)?,
    }
    eprintln!("vrsify seqmap: {} contigs", rows.len());
    Ok(())
}

#[derive(Default)]
struct VcfCounts {
    variants: u64,
    alleles: u64,
    observations: u64,
    symbolic_alts: u64,
    no_position: u64,
    unknown_contig: u64,
    non_nucleotide: u64,
    unnormalized: u64,
    unnormalized_observations: u64,
    not_fully_justified: u64,
}

fn run_convert(args: ConvertArgs) -> Result<()> {
    let vcf_path = args.vcf.as_ref().context("--vcf is required")?;
    let seqmap_path = args.seqmap.as_ref().context("--seqmap is required")?;
    let out_alleles = args.out_alleles.as_ref().context("--out-alleles is required")?;
    let out_obs = args
        .out_observations
        .as_ref()
        .context("--out-observations is required")?;
    let unidentified = Unidentified::from_args(args.variant_id_prefix.clone(), args.strict)?;

    let seqmap = load_seqmap(std::io::BufReader::new(
        File::open(seqmap_path).with_context(|| format!("opening seqmap {seqmap_path:?}"))?,
    ))?;

    let references = match &args.reference {
        Some(p) => Some(
            FastaReferences::from_path(p).with_context(|| format!("opening reference {p:?}"))?,
        ),
        None => None,
    };

    // noodles handles bgzf/plain-text detection and header parsing.
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(vcf_path)
        .with_context(|| format!("opening VCF {vcf_path:?}"))?;
    let header = reader.read_header().context("reading VCF header")?;
    // VEP writes its column layout inside the CSQ (or `vep`) INFO field's description.
    let csq_format = header
        .infos()
        .get("CSQ")
        .or_else(|| header.infos().get("vep"))
        .and_then(|m| csq_format_from_description(m.description()));
    let sample_names: Vec<String> = header.sample_names().iter().cloned().collect();

    let mut allele_out = BufWriter::new(File::create(out_alleles)?);
    let mut obs_out = BufWriter::new(File::create(out_obs)?);

    let mut seen_alleles: HashSet<String> = HashSet::new();
    let mut seen_unnormalized: HashSet<String> = HashSet::new();
    let mut c = VcfCounts::default();
    let mut missing_contigs: HashSet<String> = HashSet::new();

    for (record_no, result) in reader.records().enumerate() {
        let record = result.context("reading VCF record")?;
        let chrom = record.reference_sequence_name().to_string();
        // `chr1` vs `1`: the seqmap may be keyed either way (see `resolve_contig`), and
        // the resolved key is what the reference FASTA is then looked up by. An
        // unresolved contig is *not* a reason to drop the record — see `reject` below.
        let resolved = resolve_contig(&seqmap, &chrom);
        let pos = match record.variant_start() {
            Some(p) => usize::from(p.context("parsing POS")?) as i64,
            None => {
                // No POS at all (a telomere record, `.` in the column): there is no
                // locus to key a variant node on, so there is nothing to keep.
                c.no_position += 1;
                continue;
            }
        };
        let reference = record.reference_bases().to_string();
        let info = record.info().as_ref().to_string();

        let alts: Vec<String> = record
            .alternate_bases()
            .as_ref()
            .split(',')
            .map(|s| s.to_string())
            .filter(|s| !s.is_empty())
            .collect();

        // GT column index within FORMAT keys.
        let samples = record.samples();
        let gt_idx = samples.keys().iter().position(|k| k == "GT");

        for (i, alt) in alts.iter().enumerate() {
            let alt_number = i + 1; // VCF: 0=REF, 1=first ALT, ...
            if !is_concrete_alt(alt) {
                // Symbolic/structural ALTs (`<DEL>`, breakends) and the `*` spanning
                // deletion are an explicit scope exclusion, not an identity failure:
                // `<DEL>` is defined by INFO/END or SVLEN, which an
                // `{assembly}:{chrom}:{pos}:{ref}:{alt}` node would silently drop, and
                // `*` is not an allele of its own. Counted, never guessed at.
                c.symbolic_alts += 1;
                continue;
            }

            // Same identity policy as the MAF front end: a non-ACGTN REF or ALT (IUPAC
            // codes, caller junk) must never be given a `ga4gh:VA.` id, and neither may
            // a contig whose sequence the seqmap cannot identify.
            let reject = match resolved {
                None => {
                    missing_contigs.insert(chrom.clone());
                    c.unknown_contig += 1;
                    Some("contig is not in the seqmap".to_string())
                }
                Some(_) => None,
            }
            .or_else(|| {
                (!is_nucleotide_sequence(&reference) || !is_nucleotide_sequence(alt)).then(|| {
                    c.non_nucleotide += 1;
                    format!("non-nucleotide allele(s) '{reference}' / '{alt}'")
                })
            });

            // Per-ALT: on a multiallelic record each ALT has its own CSQ/ANN entries,
            // and a record split by `bcftools norm -m-` keeps the full original list.
            let annotation =
                extract_annotation(&info, csq_format.as_deref(), &reference, alt, alt_number);

            // Which sample carries the ALT, and from which file, is record-level
            // provenance that does not depend on the ALT earning a VRS id: the
            // observations are built the same way either way, only `variant_id` differs.
            let write_observations = |variant_id: &str,
                                          assembly: Option<String>,
                                          reference_name: &str,
                                          obs_out: &mut BufWriter<File>|
             -> Result<u64> {
                let mut n = 0;
                for (s_i, sample_name) in sample_names.iter().enumerate() {
                    let gt = match gt_idx {
                        Some(gi) => samples
                            .get_index(s_i)
                            .and_then(|s| s.as_ref().split(':').nth(gi).map(|v| v.to_string()))
                            .unwrap_or_else(|| ".".to_string()),
                        None => continue,
                    };
                    if let Some(zyg) = zygosity_for(&gt, alt_number) {
                        let obs = Observation {
                            variant_id: variant_id.to_string(),
                            sample: sample_name.clone(),
                            zygosity: zyg,
                            contig: chrom.clone(),
                            pos,
                            reference: reference.clone(),
                            alt: alt.clone(),
                            assembly: assembly.clone(),
                            reference_name: reference_name.to_string(),
                            source: args.source.clone(),
                            annotation: annotation.clone(),
                        };
                        writeln!(obs_out, "{}", serde_json::to_string(&obs.to_json())?)?;
                        n += 1;
                    }
                }
                Ok(n)
            };

            if let Some(reason) = reject {
                let Unidentified::KeepUnder(prefix) = &unidentified else {
                    bail!(
                        "VCF record {} ({chrom}:{pos} {reference}>{alt}): {reason}",
                        record_no + 1
                    );
                };
                // One assembly value for the local node identity, its metadata, and
                // every observation that references it: the seqmap when the contig
                // resolved, else the CLI's `--assembly`, else `unknown`.
                let assembly = resolved
                    .and_then(|(_, s)| s.assembly.clone())
                    .or_else(|| args.assembly.clone())
                    .unwrap_or_else(|| "unknown".to_string());
                let u = UnnormalizedVariant::new(
                    prefix, &assembly, &chrom, pos, &reference, alt, reason,
                );
                if seen_unnormalized.insert(u.key.clone()) {
                    writeln!(allele_out, "{}", serde_json::to_string(&u.to_json())?)?;
                    c.unnormalized += 1;
                }
                // The variant node is deduplicated, the observations are not: every
                // sample carrying the ALT still gets a record, pointed at that node.
                let reference_name = resolved
                    .map(|(_, s)| s.reference_name.clone())
                    .unwrap_or_else(|| chrom.clone());
                let n = write_observations(
                    &u.id(),
                    Some(assembly),
                    &reference_name,
                    &mut obs_out,
                )?;
                // Counted in the total too, so the summary matches the observation file.
                c.observations += n;
                c.unnormalized_observations += n;
                continue;
            }

            let (seqmap_key, seq) = resolved.expect("rejected above when absent");
            c.variants += 1;

            let refseq = references
                .as_ref()
                .map(|r| {
                    reference_contig(r, &[seqmap_key, chrom.as_str()])
                        .with_context(|| format!("reference FASTA is missing contig '{chrom}'"))
                })
                .transpose()?;
            let (allele, fully_justified) = match refseq {
                Some(refseq) => (
                    build_allele_normalized(seq, pos, &reference, alt, refseq)
                        .with_context(|| format!("VCF {chrom}:{pos} {reference}>{alt}"))?,
                    true,
                ),
                None => build_allele(seq, pos, &reference, alt),
            };
            if !fully_justified {
                c.not_fully_justified += 1;
            }
            let vid = allele.ga4gh_id();

            if seen_alleles.insert(vid.clone()) {
                let mut value = allele.to_output_value();
                if !fully_justified {
                    // Output-only marker; it is not part of the VRS digest.
                    value["fullyJustified"] = false.into();
                }
                writeln!(allele_out, "{}", serde_json::to_string(&value)?)?;
                c.alleles += 1;
            }

            c.observations += write_observations(
                &vid,
                seq.assembly.clone().or_else(|| args.assembly.clone()),
                &seq.reference_name,
                &mut obs_out,
            )?;
        }
    }

    allele_out.flush()?;
    obs_out.flush()?;

    eprintln!(
        "vrsify: {} variant-alleles, {} unique, {} observations, {} non-concrete ALTs skipped",
        c.variants, c.alleles, c.observations, c.symbolic_alts
    );
    if c.unnormalized > 0 {
        eprintln!(
            "vrsify: {} unnormalized variant(s) kept with {} observation(s) (unknown \
             contig {}, non-nucleotide alleles {})",
            c.unnormalized, c.unnormalized_observations, c.unknown_contig, c.non_nucleotide
        );
    }
    if c.no_position > 0 {
        eprintln!(
            "vrsify: WARNING — {} record(s) had no POS and were skipped (no locus to key \
             a variant on)",
            c.no_position
        );
    }
    if c.not_fully_justified > 0 {
        eprintln!(
            "vrsify: WARNING — {} indel allele(s) are NOT fully justified \
             (no --reference); their ga4gh ids may not match normalized variants",
            c.not_fully_justified
        );
    }
    if !missing_contigs.is_empty() {
        let mut v: Vec<_> = missing_contigs.into_iter().collect();
        v.sort();
        eprintln!(
            "vrsify: WARNING — contigs missing from seqmap (kept as unnormalized): {}",
            v.join(", ")
        );
    }
    Ok(())
}

/// Resolve a source contig name (MAF `Chromosome` / VCF `CHROM`) in a FASTA, tolerating
/// a `chr`-prefix mismatch the same way [`resolve_contig`] does for the seqmap.
fn reference_contig<'a>(
    refs: &'a FastaReferences,
    names: &[&str],
) -> Option<&'a dyn vrsify::normalize::Reference> {
    for name in names {
        for candidate in [name.to_string(), format!("chr{name}")]
            .into_iter()
            .chain(name.strip_prefix("chr").map(String::from))
        {
            if let Some(c) = refs.contig(&candidate) {
                return Some(c as &dyn vrsify::normalize::Reference);
            }
        }
    }
    None
}

/// What to do with a record that cannot be given a VRS id, settled before any work.
///
/// A `ga4gh:VA.` id is a digest: it means the same thing to everyone, so `vrsify` can
/// compute it unaided. The id of a record that *cannot* be normalized is the opposite —
/// a deterministic key over the source coordinates, unique only within whoever minted
/// it. There is no correct namespace to default to, and picking one would label a
/// stranger's data as ours, so the caller names one (or opts out with `--strict`).
///
/// The decision is made from the arguments, up front: a failed run's output files hold
/// however much was written before the error and cannot be used, so discovering a
/// missing prefix at the first unidentifiable record would throw away every record
/// converted before it. Holding the prefix in the `KeepUnder` arm also means the reject
/// path cannot reach for a namespace that was never supplied.
enum Unidentified {
    /// Keep the record as a local node, minting its id under this namespace.
    KeepUnder(String),
    /// `--strict`: fail the run instead of keeping it, so no local id is ever needed.
    Reject,
}

impl Unidentified {
    fn from_args(prefix: Option<String>, strict: bool) -> Result<Self> {
        match (strict, prefix) {
            // `--strict` keeps no local ids, so a prefix is moot rather than conflicting.
            (true, _) => Ok(Self::Reject),
            (false, Some(p)) => Ok(Self::KeepUnder(p)),
            (false, None) => bail!(
                "--variant-id-prefix is required. Records that cannot be given a VRS id \
                 (unknown contig, non-ACGTN alleles, MAF assembly mismatch) are kept \
                 under a local id, whose namespace is only meaningful to whoever minted \
                 it — there is no correct default. Pass --variant-id-prefix (e.g. \
                 'nf:variant/'), or --strict to reject such records instead of keeping \
                 them."
            ),
        }
    }
}

/// True if every base is an unambiguous/ambiguous nucleotide code VRS can carry. MAF
/// occasionally holds `NA`, VEP-style placeholders, or protein-space junk in the allele
/// columns; those rows must not be given a VRS identity.
fn is_nucleotide_sequence(bases: &str) -> bool {
    bases
        .bytes()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

#[derive(Default)]
struct MafCounts {
    rows: u64,
    variants: u64,
    alleles: u64,
    observations: u64,
    no_variant: u64,
    unknown_contig: u64,
    assembly_mismatch: u64,
    non_nucleotide: u64,
    unnormalized: u64,
    unnormalized_observations: u64,
    below_min_alt_count: u64,
    not_fully_justified: u64,
    end_position_warnings: u64,
}

fn run_maf(args: MafArgs) -> Result<()> {
    let unidentified = Unidentified::from_args(args.variant_id_prefix.clone(), args.strict)?;
    let seqmap: SeqMap = load_seqmap(std::io::BufReader::new(
        File::open(&args.seqmap)
            .with_context(|| format!("opening seqmap {:?}", args.seqmap))?,
    ))?;

    let references = match &args.reference {
        Some(p) => Some(
            FastaReferences::from_path(p).with_context(|| format!("opening reference {p:?}"))?,
        ),
        None => None,
    };

    let mut reader = std::io::BufReader::new(
        File::open(&args.maf).with_context(|| format!("opening MAF {:?}", args.maf))?,
    );
    let header = read_header(&mut reader)?;

    let mut allele_out = BufWriter::new(File::create(&args.out_alleles)?);
    let mut obs_out = BufWriter::new(File::create(&args.out_observations)?);

    let mut seen_alleles: HashSet<String> = HashSet::new();
    let mut seen_unnormalized: HashSet<String> = HashSet::new();
    let mut c = MafCounts::default();
    let mut missing_contigs: HashSet<String> = HashSet::new();
    let mut first_warnings: Vec<String> = Vec::new();

    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line).context("reading MAF row")? == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.trim().is_empty() || trimmed.starts_with('#') {
            continue;
        }
        c.rows += 1;
        let row: Vec<&str> = trimmed.split('\t').collect();

        let contig = header
            .get(&row, "Chromosome")
            .with_context(|| format!("MAF row {}: empty Chromosome", c.rows))?
            .to_string();
        let start: i64 = header
            .get(&row, "Start_Position")
            .and_then(|v| v.parse().ok())
            .with_context(|| format!("MAF row {}: unparseable Start_Position", c.rows))?;
        let end: Option<i64> = header.get(&row, "End_Position").and_then(|v| v.parse().ok());

        // MAF-form alleles, `-` placeholders preserved for the observation record.
        let ref_raw = header.get(&row, "Reference_Allele").with_context(|| {
            format!(
                "MAF row {}: missing Reference_Allele; use '-' only for an explicit empty allele",
                c.rows
            )
        })?;
        let tsa2_raw = header.get(&row, "Tumor_Seq_Allele2").with_context(|| {
            format!(
                "MAF row {}: missing Tumor_Seq_Allele2; use '-' only for an explicit empty allele",
                c.rows
            )
        })?;
        let tsa1_raw = header.get(&row, "Tumor_Seq_Allele1").unwrap_or(ref_raw);
        let reference_allele = allele_bases(Some(ref_raw));
        let Some(alt_raw) = effective_alt(ref_raw, tsa2_raw, tsa1_raw) else {
            // Both tumor alleles equal the reference: the row records no variant.
            c.no_variant += 1;
            continue;
        };
        let alt_allele = allele_bases(Some(alt_raw));

        let depth = Depth::from_row(&header, &row);
        if args.min_tumor_alt_count > 0
            && depth
                .t_alt_count
                .is_some_and(|n| n < args.min_tumor_alt_count)
        {
            c.below_min_alt_count += 1;
            continue;
        }

        let build = header
            .get(&row, "NCBI_Build")
            .map(String::from)
            .or_else(|| args.assembly.clone());

        // Assembly is part of variant *identity*: a GRCh37 row hashed against a GRCh38
        // refget accession silently mints a wrong id, so the mismatch is fatal per-row.
        let resolved = resolve_contig(&seqmap, &contig);
        let reject = match (&resolved, &build) {
            (None, _) => {
                missing_contigs.insert(contig.clone());
                c.unknown_contig += 1;
                Some("contig is not in the seqmap".to_string())
            }
            (Some((_, seq)), Some(b)) => match &seq.assembly {
                Some(sa) if !build_matches(sa, b) => {
                    c.assembly_mismatch += 1;
                    Some(format!("MAF NCBI_Build '{b}' != seqmap assembly '{sa}'"))
                }
                _ => None,
            },
            (Some(_), None) => None,
        }
        .or_else(|| {
            (!is_nucleotide_sequence(reference_allele) || !is_nucleotide_sequence(alt_allele)).then(
                || {
                    c.non_nucleotide += 1;
                    format!("non-nucleotide allele(s) '{ref_raw}' / '{alt_raw}'")
                },
            )
        });

        // Which sample, study and file a row came from is row-level provenance that does
        // not depend on the row earning a VRS id, so the observation is built the same
        // way for rejected and normalized rows; only `variant_id` differs.
        let row_no = c.rows;
        let make_obs = |variant_id: String,
                        assembly: Option<String>,
                        reference_name: String|
         -> Result<MafObservation> {
            Ok(MafObservation {
                variant_id,
                tumor_sample: header
                    .get(&row, "Tumor_Sample_Barcode")
                    .with_context(|| format!("MAF row {row_no}: empty Tumor_Sample_Barcode"))?
                    .to_string(),
                matched_normal: header
                    .get(&row, "Matched_Norm_Sample_Barcode")
                    .map(String::from),
                contig: contig.clone(),
                start,
                end,
                reference_allele: ref_raw.to_string(),
                alt_allele: alt_raw.to_string(),
                assembly,
                reference_name,
                source: args.source.clone(),
                study_id: args.study_id.clone(),
                center: header.get(&row, "Center").map(String::from),
                mutation_status: header
                    .get(&row, "Mutation_Status")
                    .map(String::from)
                    .or_else(|| Some(args.mutation_status.clone()).filter(|s| !s.is_empty())),
                annotation: MafAnnotation::from_row(&header, &row),
                depth: depth.clone(),
            })
        };

        if let Some(reason) = reject {
            let Unidentified::KeepUnder(prefix) = &unidentified else {
                bail!("MAF row {} ({contig}:{start} {ref_raw}>{alt_raw}): {reason}", c.rows);
            };
            // Use one assembly value for the local node identity, its metadata, and
            // every observation that references it. A resolved seqmap supplies the
            // assembly when the MAF and CLI do not; only an unresolved row is unknown.
            let unnormalized_assembly = build
                .clone()
                .or_else(|| resolved.and_then(|(_, s)| s.assembly.clone()))
                .unwrap_or_else(|| "unknown".to_string());
            let u = UnnormalizedVariant::new(
                prefix,
                &unnormalized_assembly,
                &contig,
                start,
                ref_raw,
                alt_raw,
                reason,
            );
            if seen_unnormalized.insert(u.key.clone()) {
                writeln!(allele_out, "{}", serde_json::to_string(&u.to_json())?)?;
                c.unnormalized += 1;
            }
            // The variant node is deduplicated, the observations are not: every sample
            // carrying the row still gets a record, pointed at the unnormalized node.
            let obs = make_obs(
                u.id(),
                Some(unnormalized_assembly),
                resolved
                    .map(|(_, s)| s.reference_name.clone())
                    .unwrap_or_else(|| contig.clone()),
            )?;
            writeln!(obs_out, "{}", serde_json::to_string(&obs.to_json())?)?;
            // Counted in the total too, so the summary matches the observation file.
            c.observations += 1;
            c.unnormalized_observations += 1;
            continue;
        }

        let (seqmap_key, seq) = resolved.expect("rejected above when absent");
        let (interval, warning) = maf_to_interbase(start, end, reference_allele, alt_allele)
            .with_context(|| format!("MAF row {}", c.rows))?;
        if let Some(w) = warning {
            c.end_position_warnings += 1;
            if first_warnings.len() < 5 {
                first_warnings.push(format!("row {}: {w}", c.rows));
            }
        }
        c.variants += 1;

        let ref_seq = references
            .as_ref()
            .map(|r| {
                reference_contig(r, &[seqmap_key, contig.as_str()]).with_context(|| {
                    format!("MAF row {}: reference FASTA is missing contig '{contig}'", c.rows)
                })
            })
            .transpose()?;

        let (allele, fully_justified) =
            build_maf_allele(seq, &interval, ref_seq).with_context(|| {
                format!("MAF row {} ({contig}:{start} {ref_raw}>{alt_raw})", c.rows)
            })?;
        if !fully_justified {
            c.not_fully_justified += 1;
        }
        let vid = allele.ga4gh_id();
        if seen_alleles.insert(vid.clone()) {
            let mut value = allele.to_output_value();
            if !fully_justified {
                // Output-only marker; it is not part of the VRS digest.
                value
                    .as_object_mut()
                    .expect("allele serializes to an object")
                    .insert("fullyJustified".into(), false.into());
            }
            writeln!(allele_out, "{}", serde_json::to_string(&value)?)?;
            c.alleles += 1;
        }

        let obs = make_obs(
            vid,
            build.or_else(|| seq.assembly.clone()),
            seq.reference_name.clone(),
        )?;
        writeln!(obs_out, "{}", serde_json::to_string(&obs.to_json())?)?;
        c.observations += 1;
    }

    allele_out.flush()?;
    obs_out.flush()?;

    eprintln!(
        "vrsify maf: {} rows → {} variant-alleles, {} unique alleles, {} observations",
        c.rows, c.variants, c.alleles, c.observations
    );
    if c.below_min_alt_count > 0 {
        eprintln!(
            "vrsify maf: {} rows dropped for t_alt_count < {}",
            c.below_min_alt_count, args.min_tumor_alt_count
        );
    }
    if c.no_variant > 0 {
        eprintln!("vrsify maf: {} rows had no non-reference tumor allele", c.no_variant);
    }
    if c.unnormalized > 0 {
        eprintln!(
            "vrsify maf: {} unnormalized variant(s) kept with {} observation(s) (unknown \
             contig {}, assembly mismatch {}, non-nucleotide alleles {})",
            c.unnormalized,
            c.unnormalized_observations,
            c.unknown_contig,
            c.assembly_mismatch,
            c.non_nucleotide
        );
    }
    if !missing_contigs.is_empty() {
        let mut v: Vec<_> = missing_contigs.into_iter().collect();
        v.sort();
        eprintln!("vrsify maf: contigs missing from seqmap: {}", v.join(", "));
    }
    if c.not_fully_justified > 0 {
        eprintln!(
            "vrsify maf: WARNING — {} indel allele(s) are NOT fully justified (no \
             --reference); their ga4gh ids will not match vrs-python/ClinVar",
            c.not_fully_justified
        );
    }
    if c.end_position_warnings > 0 {
        eprintln!(
            "vrsify maf: {} row(s) had an End_Position inconsistent with \
             Start_Position + Reference_Allele (Reference_Allele used):",
            c.end_position_warnings
        );
        for w in &first_warnings {
            eprintln!("  {w}");
        }
    }
    Ok(())
}
