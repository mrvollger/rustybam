use clap::{load_yaml, App, AppSettings};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustybam::bamstats;
use rustybam::bed;
use rustybam::nucfreq;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let app = App::from(yaml).setting(AppSettings::SubcommandRequiredElseHelp);
    let matches = app.get_matches();

    if let Some(matches) = matches.subcommand_matches("stats") {
        run_stats(matches);
    } else if let Some(matches) = matches.subcommand_matches("nucfreq") {
        run_nucfreq(matches);
    }
}

pub fn run_stats(args: &clap::ArgMatches) {
    // parse arguments
    let threads = args.value_of_t("threads").unwrap_or(8);
    eprintln!("Number of threads: {}", threads);
    let mut bam = match args.value_of("BAM") {
        Some(bam_f) => {
            bam::Reader::from_path(bam_f).unwrap_or_else(|_| panic!("Failed to open {}", bam_f))
        }
        _ => bam::Reader::from_stdin().unwrap(),
    };
    let qbed = args.is_present("qbed");

    // open bam
    bam.set_threads(threads).unwrap();
    let bam_header = bam::Header::from_template(bam.header());

    // get stats
    bamstats::print_cigar_stats_header(qbed);
    for (idx, rec) in bam.records().enumerate() {
        eprint!("\rProccesing: {}", idx);
        let stats = bamstats::cigar_stats(rec.unwrap());
        bamstats::print_cigar_stats(stats, qbed, &bam_header);
    }
    eprintln!();
}

pub fn run_nucfreq(args: &clap::ArgMatches) {
    // parse arguments
    let threads = args.value_of_t("threads").unwrap_or(8);
    eprintln!("Number of threads: {}", threads);
    let bam_f = args.value_of("BAM").unwrap();
    let mut bam =
        bam::IndexedReader::from_path(bam_f).unwrap_or_else(|_| panic!("Failed to open {}", bam_f));

    nucfreq::print_nucfreq_header();
    // nuc freq on region
    if args.is_present("region") {
        let rgn = bed::parse_region(args.value_of("region").unwrap());
        let vec = nucfreq::region_nucfreq(&mut bam, &rgn);
        nucfreq::print_nucfreq(vec, &rgn);
    }
    //nucfreq on bed
    if args.is_present("bed") {
        let bed_f = args.value_of("bed").expect("Unable to read bedfile");
        for rgn in bed::parse_bed(bed_f) {
            let vec = nucfreq::region_nucfreq(&mut bam, &rgn);
            nucfreq::print_nucfreq(vec, &rgn);
        }
    }
}
