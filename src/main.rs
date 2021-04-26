use clap::{load_yaml, App, AppSettings};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustybam::bed;
use rustybam::nucfreq;
use rustybam::paf;
use rustybam::suns;
use rustybam::{bamstats, paf::trim_paf_to_rgn};
use std::fs;
use std::io::{self, BufRead, BufReader};
use std::time::Instant;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let app = App::from(yaml).setting(AppSettings::SubcommandRequiredElseHelp);
    let matches = app.get_matches();

    if let Some(matches) = matches.subcommand_matches("stats") {
        run_stats(matches);
    } else if let Some(matches) = matches.subcommand_matches("nucfreq") {
        run_nucfreq(matches);
    } else if let Some(matches) = matches.subcommand_matches("suns") {
        run_suns(matches);
    } else if let Some(matches) = matches.subcommand_matches("liftover") {
        run_liftover(matches);
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
    // set the number of threads
    let threads = args.value_of_t("threads").unwrap_or(8);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    eprintln!("Number of threads: {}", threads);

    // read the bam
    let bam_f = args
        .value_of("BAM")
        .expect("Must provide an indexed alignment file (bam/cram)");

    // add the nuc freq regions to process.
    let mut rgns = Vec::new();
    if args.is_present("region") {
        rgns.push(bed::parse_region(args.value_of("region").unwrap()));
    }
    if args.is_present("bed") {
        let bed_f = args.value_of("bed").expect("Unable to read bedfile");
        rgns.append(&mut bed::parse_bed(bed_f));
    }

    for rgn in rgns {
        // say the max window size a region can be before printing
        let med_rgns = bed::split_region(&rgn, 1_000_000);
        // split the windows into windows of that size
        for med_rgn in med_rgns {
            let small_rgns = bed::split_region(&med_rgn, 10_000);
            // generate the nucfreqs
            let vec: Vec<nucfreq::Nucfreq> = small_rgns
                .into_par_iter()
                .map(|r| nucfreq::region_nucfreq(bam_f, &r, 4))
                .flatten()
                .collect();

            // print the results
            if args.is_present("small") {
                nucfreq::small_nucfreq(&vec)
            } else {
                nucfreq::print_nucfreq_header();
                nucfreq::print_nucfreq(&vec);
            }
        }
    }
}

pub fn run_suns(args: &clap::ArgMatches) {
    let kmer_size = args.value_of_t("kmersize").unwrap_or(21);
    let max_interval = args.value_of_t("maxsize").unwrap_or(std::usize::MAX);
    let fastafile = args.value_of("fasta").expect("Fasta file reuqired!");
    let genome = suns::Genome::from_file(fastafile);
    let sun_intervals = genome.find_sun_intervals(kmer_size);
    println!("#chr\tstart\tend\tsun_seq");
    for (chr, start, end, seq) in &sun_intervals {
        if end - start < max_interval {
            println!(
                "{}\t{}\t{}\t{}",
                chr,
                start,
                end,
                std::str::from_utf8(seq).unwrap()
            );
        }
    }
    if args.is_present("validate") {
        suns::validate_suns(&genome, &sun_intervals, kmer_size);
    }
}

pub fn run_liftover(args: &clap::ArgMatches) {
    let paf_buff: Box<dyn BufRead> = match args.value_of("paf") {
        Some(paf_f) => Box::new(BufReader::new(fs::File::open(paf_f).unwrap())),
        _ => Box::new(BufReader::new(io::stdin())),
    };
    let start = Instant::now();
    let bed = args.value_of("bed").expect("Bed file required!");
    let paf = paf::read_paf(paf_buff);
    let rgns = bed::parse_bed(bed);
    let duration = start.elapsed();
    eprintln!("Time elapsed reading paf and bed: {:.3?}", duration);
    // whether the input bed is for the query.
    let mut invert_query = false;
    if args.is_present("querybed") {
        invert_query = true;
    }
    //
    let start = Instant::now();
    for mut rgn in rgns {
        let new_paf = trim_paf_to_rgn(&rgn, &paf, invert_query);
        for rec in new_paf {
            if rgn.id == "None" {
                rgn.id = format!("{}_{}_{}", rgn.name, rgn.st, rgn.en)
            }
            println!("{}\tid:Z:{}", rec, rgn.id);
        }
    }
    let duration = start.elapsed();
    eprintln!("Time elapsed during liftover: {:.3?}", duration);
}
