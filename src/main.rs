use clap::{App, load_yaml};
use rustybam::bamstats;
use rustybam::nucfreq;
use rust_htslib::bam;
use rust_htslib::bam::Read;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from(yaml).get_matches();

    if let Some(matches) = matches.subcommand_matches("stats"){
        run_stats(matches);
    } else if let Some(matches) = matches.subcommand_matches("nucfreq"){
        run_nucfreq(matches);
    }
}


pub fn run_stats( args : &clap::ArgMatches){
    // parse arguments 
    let threads = args.value_of_t("threads").unwrap_or(8);
    eprint!("Number of threads: {}\n", threads);
    let mut bam = match args.value_of("BAM") {
        Some(bam_f) => bam::Reader::from_path(bam_f).expect(&format!("Failed to open {}.", bam_f)),
        _ => bam::Reader::from_stdin().unwrap()
    };
    let qbed = args.is_present("qbed");

    // open bam
    bam.set_threads(threads).unwrap();
    let bam_header = bam::Header::from_template(bam.header());

    // get stats 
    bamstats::print_cigar_stats_header(qbed);
    for (idx, rec) in bam.records().enumerate(){
        eprint!("\rProccesing: {}", idx);
        let stats = bamstats::cigar_stats(rec.unwrap());
        bamstats::print_cigar_stats(stats, qbed, &bam_header);
   }
   eprint!("\n");
}

pub fn run_nucfreq( args : &clap::ArgMatches){
    // parse arguments 
    let threads = args.value_of_t("threads").unwrap_or(8);
    eprint!("Number of threads: {}\n", threads);
    let bam_f = args.value_of("BAM").unwrap();
    let region = args.value_of("region").unwrap();
    let mut bam = bam::IndexedReader::from_path(bam_f).expect(&format!("Failed to open {}.", bam_f));

    let vec = nucfreq::region_nucfreq(&mut bam, region);
    println!("{}", vec.len());
    nucfreq::print_nucfreq(vec);
}

