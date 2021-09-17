use bio::io::fasta;
use bio::io::fastq;
use clap::{crate_version, load_yaml, App, AppSettings};
use itertools::Itertools;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustybam::bamstats;
use rustybam::bed;
use rustybam::liftover;
use rustybam::myio;
use rustybam::nucfreq;
use rustybam::paf;
use rustybam::suns;
use std::io;
use std::time::Instant;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let app = App::from(yaml)
        .version(crate_version!())
        .setting(AppSettings::SubcommandRequiredElseHelp);
    let matches = app.get_matches();

    let threads = matches.value_of_t("threads").unwrap_or(8);
    std::env::set_var("RAYON_NUM_THREADS", threads.to_string());

    if let Some(matches) = matches.subcommand_matches("stats") {
        run_stats(matches);
    } else if let Some(matches) = matches.subcommand_matches("nucfreq") {
        run_nucfreq(matches);
    } else if let Some(matches) = matches.subcommand_matches("suns") {
        run_suns(matches);
    } else if let Some(matches) = matches.subcommand_matches("liftover") {
        run_liftover(matches);
    } else if let Some(matches) = matches.subcommand_matches("repeat") {
        run_longest_repeats(matches);
    } else if let Some(matches) = matches.subcommand_matches("bedlength") {
        run_bedlength(matches);
    } else if let Some(matches) = matches.subcommand_matches("breakpaf") {
        run_break_paf(matches);
    } else if let Some(matches) = matches.subcommand_matches("fastq-split") {
        run_split_fastq(matches);
    } else if let Some(matches) = matches.subcommand_matches("fasta-split") {
        run_split_fasta(matches);
    }
}

pub fn run_stats(args: &clap::ArgMatches) {
    // parse arguments
    let threads = args.value_of_t("threads").unwrap_or(8);
    eprintln!("Number of threads: {}", threads);
    let qbed = args.is_present("qbed");
    let paf = args.is_present("paf");
    bamstats::print_cigar_stats_header(qbed);

    if paf {
        let file = args.value_of("BAM").unwrap_or("-");
        let mut idx = 1;
        for paf in paf::Paf::from_file(file).records {
            eprint!("\rProcessing: {}", idx);
            let stats = bamstats::stats_from_paf(paf);
            bamstats::print_cigar_stats(stats, qbed);
            idx += 1;
        }
        eprintln!();
        return;
    }
    // Not a paf so lets read in the bam

    // we want to do bam reading
    let mut bam = match args.value_of("BAM") {
        Some(bam_f) => {
            bam::Reader::from_path(bam_f).unwrap_or_else(|_| panic!("Failed to open {}", bam_f))
        }
        _ => bam::Reader::from_stdin().unwrap(),
    };

    // open bam
    bam.set_threads(threads).unwrap();
    let bam_header = bam::Header::from_template(bam.header());

    // get stats
    for (idx, rec) in bam.records().enumerate() {
        eprint!("\rProcessing: {}", idx + 1);
        let stats = bamstats::cigar_stats(rec.unwrap(), &bam_header);
        bamstats::print_cigar_stats(stats, qbed);
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
    let fastafile = args.value_of("fasta").expect("Fasta file required!");
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

pub fn run_longest_repeats(args: &clap::ArgMatches) {
    let minsize = args.value_of_t("min").unwrap_or(21);
    let fastafile = args.value_of("fasta").expect("Fasta file required!");
    let genome = suns::Genome::from_file(fastafile);
    let unique_intervals = genome.get_longest_perfect_repeats(minsize);
    println!("#chr\tstart\tend\trepeat_length");
    for (chr, start, length) in &unique_intervals {
        println!("{}\t{}\t{}\t{}", chr, start, start + length, length - 1,);
    }
}

pub fn run_liftover(args: &clap::ArgMatches) {
    let start = Instant::now();
    // read in the bed
    let bed = args.value_of("bed").expect("Bed file required!");
    let rgns = bed::parse_bed(bed);
    // read in the file
    let paf_file = args.value_of("paf").unwrap_or("-");
    let paf = paf::Paf::from_file(paf_file);
    // end timer
    let duration = start.elapsed();
    eprintln!("Time elapsed reading paf and bed: {:.3?}", duration);
    // whether the input bed is for the query.
    let mut invert_query = false;
    if args.is_present("qbed") {
        invert_query = true;
    }
    //
    let start = Instant::now();
    let new_recs = liftover::trim_paf_by_rgns(&rgns, &paf.records, invert_query);

    // if largest set report only the largest alignment for the record
    let largest = args.is_present("largest");
    if largest {
        for (_key, group) in &new_recs
            .into_iter()
            .sorted_by_key(|pac_rec| pac_rec.id.clone())
            .group_by(|paf_rec| paf_rec.id.clone())
        {
            let largest_rec = group.max_by_key(|p| (p.t_en - p.t_st)).unwrap();
            println!("{}", largest_rec);
        }
    } else {
        for rec in new_recs {
            println!("{}", rec);
        }
    }

    let duration = start.elapsed();
    eprintln!("Time elapsed during liftover: {:.3?}", duration);
}

pub fn run_bedlength(args: &clap::ArgMatches) {
    let start = Instant::now();
    // read in the bed
    let bed = args.value_of("bed").expect("Bed file required!");
    let rgns = bed::parse_bed(bed);
    let count: u64 = rgns.into_iter().map(|rgn| rgn.en - rgn.st).sum();
    println!("{}", count);
    let duration = start.elapsed();
    eprintln!("Time elapsed during bedlength: {:.3?}", duration);
}

pub fn run_break_paf(args: &clap::ArgMatches) {
    let start = Instant::now();
    // break length
    let break_length = args.value_of_t("max-size").unwrap_or(100);
    // read in the file
    let paf_file = args.value_of("paf").unwrap_or("-");
    let paf = paf::Paf::from_file(paf_file);

    for mut paf in paf.records {
        paf.aligned_pairs();
        let pafs = liftover::break_paf_on_indels(&paf, break_length);
        for trimed_paf in pafs {
            println!("{}", trimed_paf);
        }
    }
    // end timer
    let duration = start.elapsed();
    eprintln!("Time elapsed breaking paf on indels: {:.3?}", duration);
}

pub fn run_split_fastq(args: &clap::ArgMatches) {
    let start = Instant::now();
    let files: Vec<_> = args.values_of("fastq").unwrap().collect();

    let mut outs = Vec::new();
    for f in files {
        let handle = myio::writer(f);
        outs.push(fastq::Writer::new(handle));
    }

    let mut records = fastq::Reader::new(io::stdin()).records();
    let mut out_idx = 0;
    while let Some(Ok(record)) = records.next() {
        outs[out_idx]
            .write_record(&record)
            .expect("Error writing record.");
        out_idx += 1;
        if out_idx == outs.len() {
            out_idx = 0;
        }
    }
    // end timer
    let duration = start.elapsed();
    eprintln!("Time elapsed splitting fastq: {:.3?}", duration);
}

pub fn run_split_fasta(args: &clap::ArgMatches) {
    let start = Instant::now();
    let files: Vec<_> = args.values_of("fasta").unwrap().collect();

    let mut outs = Vec::new();
    for f in files {
        let handle = myio::writer(f);
        outs.push(fasta::Writer::new(handle));
    }

    let mut records = fasta::Reader::new(io::stdin()).records();
    let mut out_idx = 0;
    while let Some(Ok(record)) = records.next() {
        outs[out_idx]
            .write_record(&record)
            .expect("Error writing record.");
        out_idx += 1;
        if out_idx == outs.len() {
            out_idx = 0;
        }
    }
    // end timer
    let duration = start.elapsed();
    eprintln!("Time elapsed splitting fasta: {:.3?}", duration);
}
