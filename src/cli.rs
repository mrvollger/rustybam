use clap::IntoApp;
use clap::{App, AppSettings, Parser, Subcommand};

#[derive(Parser, Debug)]
#[clap(global_setting(AppSettings::PropagateVersion))]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(global_setting(AppSettings::InferSubcommands))]
#[clap(global_setting(AppSettings::HelpExpected))]
#[clap(global_setting(AppSettings::UseLongFormatForHelpSubcommand))]
#[clap(setting(AppSettings::SubcommandRequiredElseHelp))]
#[clap(author, version, about)]
pub struct Cli {
    /// Threads for decompression.
    #[clap(short, long, default_value_t = 8)]
    pub threads: usize,

    /// Logging level [-v: Info, -vv: Debug, -vvv: Trace].
    #[clap(short, long, parse(from_occurrences), help_heading = "DEBUG")]
    pub verbose: usize,

    #[clap(subcommand)]
    pub command: Option<Commands>,
}

///
/// This structure contains all the subcommands for rustybam and their help descriptions.
///
/// Because of naming conventions for rust enums the commands names have
/// different capitalization than on the command line.
/// For example, the `Liftover` enum is invoked using `rustybam liftover`
/// and the `TrimPaf` command with `rustybam trim-paf`.
///
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Get percent identity stats from a sam/bam/cram or PAF.
    ///
    /// Requires =/X operations in the CIGAR string!
    Stats {
        /// sam/bam/cram/file
        #[clap(default_value = "-")]
        bam: String,
        /// print query coordinates first
        #[clap(short, long)]
        qbed: bool,
        /// The input is paf format
        /// (must have cg tag with extended cigar or cs tag).
        #[clap(short, long)]
        paf: bool,
    },
    /// Count the number of bases in a bed file.
    #[clap(visible_aliases = &["bedlen", "bl", "bedlength"])]
    BedLength {
        /// a bed file
        #[clap(default_value = "-")]
        bed: String,
        /// make the output human readable (Mbp)
        #[clap(short, long)]
        readable: bool,
        /// count bases for each category in this column <COLUMN>.
        #[clap(short, long)]
        column: Option<u8>,
    },
    /// Filter PAF records in various ways.
    Filter {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// minimum number of aligned bases across all alignments between a target and query in order to keep those records
        #[clap(short, long, default_value_t = 0)]
        paired_len: u64,
        /// minimum alignment length
        #[clap(short, long, default_value_t = 0)]
        aln: u64,
        /// minimum query length
        #[clap(short, long, default_value_t = 0)]
        query: u64,
    },
    /// Invert the target and query sequences in a PAF along with the CIGAR string.
    Invert {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
    },
    /// Liftover target sequence coordinates onto query sequence using a PAF.
    ///
    /// This is a function for lifting over coordinates from a reference (<BED>) to a query using a PAF file from minimap2 or unimap (note, you can use `paftools.js sam2paf` to convert SAM data to PAF format).
    /// The returned file is a PAF file that is trimmed to the regions in the bed file. Even the cigar in the returned PAF file is trimmed so it can be used downstream! Additionally, a tag with the format `id:Z:<>` is added to the PAF where `<>` is either the 4th column of the input bed file or if not present `chr_start_end`.
    Liftover {
        /// PAF file from minimap2 or unimap run with -c and --eqx [i.e. the PAF file must have the cg tag and use extended CIGAR opts (=/X)].
        #[clap(default_value = "-")]
        paf: String,
        /// BED file of reference regions to liftover to the query.
        #[clap(short, long)]
        bed: String,
        /// Specifies that the BED file contains query coordinates to be lifted onto the reference (reverses direction of liftover).
        ///
        /// Note, that this will make the query in the input `PAF` the target in the output `PAF`.
        #[clap(short, long)]
        qbed: bool,
        /// If multiple records overlap the same region in the <bed> return only the largest liftover. The default is to return all liftovers.
        #[clap(short, long)]
        largest: bool,
    },
    /// Trim paf records that overlap in query sequence.
    ///
    /// This idea is to mimic some of the trimming that happens in PAV to improve breakpoint detection. Starts with the largest overlap and iterates until no query overlaps remain.
    #[clap(visible_aliases = &["trim", "tp"])]
    TrimPaf {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// Value added for a matching base.
        #[clap(short, long, default_value_t = 1)]
        match_score: i32,
        /// Value subtracted for a mismatching base.
        #[clap(short, long, default_value_t = 1)]
        diff_score: i32,
        /// Value subtracted for an insertion or deletion.
        #[clap(short, long, default_value_t = 1)]
        indel_score: i32,
    },
    /// Orient paf records so that most of the bases are in the forward direction.
    ///
    /// Optionally scaffold the queriers so that there is one query per target.
    Orient {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// Generate ~1 query per target that scaffolds together all the records that map to that target sequence.
        ///
        /// The order of the scaffold will be determined by the average target position across all the queries, and the name of the new scaffold will be
        #[clap(short, long)]
        scaffold: bool,
        /// space to add between records
        #[clap(short, long, default_value_t = 1_000_000)]
        insert: u64,
    },
    /// Break PAF records with large indels into multiple records (useful for SafFire).
    #[clap(visible_aliases = &["breakpaf", "bp"])]
    BreakPaf {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// maximum indel size to keep in the paf
        #[clap(short, long, default_value_t = 100)]
        max_size: u32,
    },
    /// Convert a PAF file into a SAM file. Warning, all alignments will be marked as primary!
    #[clap(visible_aliases = &["paftosam", "p2s", "paf2sam"])]
    PafToSam {
        /// PAF file from minimap2 or unimap.
        #[clap(default_value = "-")]
        paf: String,
    },
    /// Reads in a fastx from stdin and divides into files (can compress by adding .gz).
    /// Input can be either compressed or uncompressed!
    #[clap(visible_aliases = &["fxsplit", "fasta-split", "fastq-split" ,"fa-split", "fq-split"])]
    FastxSplit {
        /// list of fasta files
        fasta: Vec<String>,
    },
    /// Mimic bedtools getfasta but allow for bgzip in both bed and fasta inputs.
    #[clap(visible_aliases = &["getfasta", "gf"])]
    GetFasta {
        /// fasta file to extract sequences from
        #[clap(short, long, default_value = "-")]
        fasta: String,
        /// bed file of regions to extract
        #[clap(short, long)]
        bed: String,
        /// revcomp sequence if the strand is "-"
        #[clap(short, long)]
        strand: bool,
        /// add the name (4th column) to the header of the fasta output
        #[clap(short, long)]
        name: bool,
    },
    /// Get the frequencies of each bp at each position.
    Nucfreq {
        /// sam/bam/cram/file
        #[clap(default_value = "-")]
        bam: String,
        /// print nucfreq info from the input region e.g "chr1:1-1000"
        #[clap(short, long)]
        region: Option<String>,
        /// print nucfreq info from regions in the bed file
        /// output is optionally tagged using the 4th column
        #[clap(short, long)]
        bed: Option<String>,
        /// smaller output format
        #[clap(short, long)]
        small: bool,
    },
    /// Report the longest exact repeat length at every position in a fasta.
    Repeat {
        /// a fasta file
        #[clap(default_value = "-")]
        fasta: String,
        /// The smallest repeat length to report
        #[clap(short, long, default_value_t = 21)]
        min: usize,
    },
    /// Extract the intervals in a genome (fasta) that are made up of SUNs.
    Suns {
        /// fasta file with the genome
        #[clap(short, long, default_value = "-")]
        fasta: String,
        /// the size of the required unique kmer
        #[clap(short, long, default_value_t = 21)]
        kmer_size: usize,
        /// the maximum size SUN interval to report
        #[clap(short, long, default_value_t = std::usize::MAX)]
        max_size: usize,
        /// confirm all the SUNs (very slow) only for debugging.
        #[clap(short, long)]
        validate: bool,
    },
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> App<'static> {
    Cli::into_app()
}
