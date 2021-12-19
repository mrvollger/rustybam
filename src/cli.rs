use clap::IntoApp;
use clap::{App, AppSettings, Parser, Subcommand};

#[derive(Parser, Debug)]
#[clap(author, version, about)]
#[clap(global_setting(AppSettings::PropagateVersion))]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(global_setting(AppSettings::UseLongFormatForHelpSubcommand))]
#[clap(setting(AppSettings::SubcommandRequiredElseHelp))]
pub struct Cli {
    /// threads for decompression
    #[clap(short, long, default_value_t = 8)]
    pub threads: usize,

    /// logging level.
    #[clap(short, long, parse(from_occurrences))]
    pub verbose: usize,

    #[clap(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// get percent identity stats from a sam/bam/cram or PAF
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

    /// count bases in a bed file
    #[clap(visible_aliases = &["bedlen", "bl", "bedlength"])]
    BedLength {
        /// a bed file
        #[clap(default_value = "-")]
        bed: String,
        /// make the output human readable (Mbp)
        #[clap(short, long)]
        readable: bool,
    },
    /// filter PAF records in various ways
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
    ///  invert the target and query sequences in a PAF along with the cg tag.
    Invert {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
    },
    ///  liftover target sequence coordinates onto query sequence using a PAF
    Liftover {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// bed file of regions to liftover
        #[clap(short, long)]
        bed: String,
        /// the bed contains query coordinates to be lifted (note the query in the original PAF will become the target in the output)
        #[clap(short, long)]
        qbed: bool,
        /// the bed contains query coordinates to be lifted (note the query in the original PAF will become the target in the output)
        #[clap(short, long)]
        largest: bool,
    },
    /// orient paf records so that most of the bases are in the forward direction
    Orient {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// make fake query names that scaffold together all the records that map to one target sequence.
        /// the order of the scaffold will be determined by the middle position of the largest alignment.
        #[clap(short, long)]
        scaffold: bool,
        /// space to add between records
        #[clap(short, long, default_value_t = 1_000_000)]
        insert: u64,
    },
    /// break up paf on indels of a certain size
    #[clap(visible_aliases = &["breakpaf", "bp"])]
    BreakPaf {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// maximum indel size to keep in the paf
        #[clap(short, long, default_value_t = 100)]
        max_size: u32,
    },
    /// reads in a fasta from stdin and divides into files (can compress by adding .gz)
    #[clap(visible_aliases = &["fastasplit", "fasplit"])]
    FastaSplit {
        /// list of fasta files
        fasta: Vec<String>,
    },
    /// reads in a fastq from stdin and divides into files (can compress by adding .gz)
    #[clap(visible_aliases = &["fastqsplit", "fqsplit"])]
    FastqSplit {
        /// list of fastq files
        fastq: Vec<String>,
    },
    /// mimic bedtools getfasta but allow for bgzip in both bed and fasta inputs.
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
    /// get the frequencies of each bp at each position.
    Nucfreq {
        // sam/bam/cram/file
        #[clap(default_value = "-")]
        bam: String,
        /// print nucfreq info from the input region e.g "chr1:1-1000"
        #[clap(short, long)]
        region: Option<String>,
        /// print nucfreq info from regions in the bed file
        /// output is optionally tagged using the 4th column
        #[clap(short, long)]
        bed: Option<String>,
        // smaller output format
        #[clap(short, long)]
        small: bool,
    },
    /// report the longest repeat length at every position in a fasta
    Repeat {
        // a fasta file
        #[clap(default_value = "-")]
        fasta: String,
        /// The smallest repeat length to report
        #[clap(short, long, default_value_t = 21)]
        min: usize,
    },
    /// extract the intervals in a genome (fasta) that are made up of SUNs
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
