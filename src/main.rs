use colored::Colorize;
use env_logger::{Builder, Target};
use itertools::Itertools;
use log::LevelFilter;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::faidx;
use rustybam::cli::Commands;
use rustybam::fastx;
use rustybam::paf::paf_swap_query_and_target;
use rustybam::*;
use std::collections::HashMap;
use std::time::Instant;

fn main() {
    parse_cli();
}

pub fn parse_cli() {
    let pg_start = Instant::now();
    let args = cli::make_cli_parse();
    let matches = cli::make_cli_app().get_matches();
    let subcommand = matches.subcommand_name().unwrap();

    // set the logging level
    let min_log_level = match matches.occurrences_of("verbose") {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    Builder::new()
        .target(Target::Stderr)
        .filter(None, min_log_level)
        .init();

    log::debug!("DEBUG logging enabled");
    log::trace!("TRACE logging enabled");

    // set up number of threads to use globally
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    match &args.command {
        //
        // Run Stats
        //
        Some(Commands::Stats { bam, qbed, paf }) => {
            bamstats::print_cigar_stats_header(*qbed);
            if *paf {
                for paf in paf::Paf::from_file(bam).records {
                    let stats = bamstats::stats_from_paf(paf);
                    bamstats::print_cigar_stats(stats, *qbed);
                }
                return;
            }
            // Not a paf so lets read in the bam
            let mut bam_reader = if bam == "-" {
                bam::Reader::from_stdin().unwrap()
            } else {
                bam::Reader::from_path(bam).unwrap_or_else(|_| panic!("Failed to open {}", bam))
            };

            // open bam
            bam_reader.set_threads(args.threads).unwrap();
            let bam_header = bam::Header::from_template(bam_reader.header());

            // get stats
            for (_idx, rec) in bam_reader.records().enumerate() {
                let rec = rec.unwrap();
                if !rec.is_unmapped() {
                    let stats = bamstats::cigar_stats(rec, &bam_header);
                    bamstats::print_cigar_stats(stats, *qbed);
                }
            }
        }
        //
        // Run Nucfreq
        //
        Some(Commands::Nucfreq {
            bam,
            region,
            bed,
            small,
        }) => {
            // add the nuc freq regions to process.
            let mut rgns = Vec::new();
            // add regions
            if let Some(region_f) = region {
                rgns.push(bed::parse_region(region_f));
            }
            // add bed
            if let Some(bed_f) = bed {
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
                        .map(|r| nucfreq::region_nucfreq(bam, &r, 4))
                        .flatten()
                        .collect();

                    // print the results
                    if *small {
                        nucfreq::small_nucfreq(&vec)
                    } else {
                        nucfreq::print_nucfreq_header();
                        nucfreq::print_nucfreq(&vec);
                    }
                }
            }
        }
        //
        // Run Repeat
        //
        Some(Commands::Repeat { fasta, min }) => {
            let genome = suns::Genome::from_file(fasta);
            let unique_intervals = genome.get_longest_perfect_repeats(*min);
            println!("#chr\tstart\tend\trepeat_length");
            for (chr, start, length) in &unique_intervals {
                println!("{}\t{}\t{}\t{}", chr, start, start + length, length - 1,);
            }
        }
        //
        // Run Suns
        //
        Some(Commands::Suns {
            fasta,
            kmer_size,
            max_size,
            validate,
        }) => {
            let genome = suns::Genome::from_file(fasta);
            let sun_intervals = genome.find_sun_intervals(*kmer_size);
            println!("#chr\tstart\tend\tsun_seq");
            for (chr, start, end, seq) in &sun_intervals {
                if end - start < *max_size {
                    println!(
                        "{}\t{}\t{}\t{}",
                        chr,
                        start,
                        end,
                        std::str::from_utf8(seq).unwrap()
                    );
                }
            }
            if *validate {
                suns::validate_suns(&genome, &sun_intervals, *kmer_size);
            }
        }
        //
        // Run Bedlength
        //
        Some(Commands::BedLength {
            bed,
            readable,
            column,
        }) => {
            let rgns = bed::parse_bed(bed);
            match *column {
                Some(c) => {
                    let mut dict = HashMap::new();
                    for rgn in rgns {
                        let bp = dict.entry(rgn.get_column(c).clone()).or_insert(0_u64);
                        *bp += rgn.en - rgn.st;
                    }
                    log::trace!("{:?}", dict);
                    for (key, count) in dict {
                        if *readable {
                            println!("{}\t{}", key, (count as f64) / 1e6);
                        } else {
                            println!("{}\t{}", key, count);
                        }
                    }
                }
                None => {
                    let count: u64 = rgns.into_iter().map(|rgn| rgn.en - rgn.st).sum();
                    if *readable {
                        println!("{}", (count as f64) / 1e6);
                    } else {
                        println!("{}", count);
                    }
                }
            }
        }
        //
        // Run Invert
        //
        Some(Commands::Invert { paf }) => {
            let paf = paf::Paf::from_file(paf);
            for rec in &paf.records {
                let inv_rec = paf_swap_query_and_target(rec);
                println!("{}", inv_rec);
            }
        }
        //
        // Run Liftover
        //
        Some(Commands::Liftover {
            paf,
            bed,
            qbed,
            largest,
        }) => {
            let rgns = bed::parse_bed(bed);
            // read in the file
            let paf = paf::Paf::from_file(paf);

            // trim the records
            let new_recs = liftover::trim_paf_by_rgns(&rgns, &paf.records, *qbed);

            // if largest set report only the largest alignment for the record
            if *largest {
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
        }
        //
        // Run TrimPaf
        //
        Some(Commands::TrimPaf {
            paf,
            match_score,
            diff_score,
            indel_score,
            remove_contained,
        }) => {
            let mut paf = paf::Paf::from_file(paf);
            paf.overlapping_paf_recs(*match_score, *diff_score, *indel_score, *remove_contained);
            for rec in &paf.records {
                println!("{}", rec);
            }
        }
        //
        // Run Filter
        //
        Some(Commands::Filter {
            paf,
            paired_len,
            aln,
            query,
        }) => {
            let mut paf = paf::Paf::from_file(paf);
            log::info!("{} PAF records BEFORE filtering.", paf.records.len());
            paf.filter_query_len(*query);
            paf.filter_aln_len(*aln);
            paf.filter_aln_pairs(*paired_len);
            log::info!("{} PAF records AFTER filtering.", paf.records.len());
            for rec in paf.records {
                println!("{}", rec);
            }
        }
        //
        // Run Orient
        //
        Some(Commands::Orient {
            paf,
            scaffold,
            insert,
        }) => {
            //orient_records(paf, *scaffold, *insert);
            let mut paf = paf::Paf::from_file(paf);
            paf.orient();
            if *scaffold {
                paf.scaffold(*insert);
            }
            for rec in &paf.records {
                println!("{}", rec);
            }
        }
        //
        // Run Breakpaf
        //
        Some(Commands::BreakPaf { paf, max_size }) => {
            // read in the file
            let paf = paf::Paf::from_file(paf);
            for mut paf in paf.records {
                paf.aligned_pairs();
                let pafs = liftover::break_paf_on_indels(&paf, *max_size);
                for trimed_paf in pafs {
                    println!("{}", trimed_paf);
                }
            }
        }
        //
        // Run PafToSam
        //
        Some(Commands::PafToSam { paf, fasta }) => {
            let fasta_reader = fasta
                .as_ref()
                .map(|fasta| faidx::Reader::from_path(fasta).unwrap());

            let paf = paf::Paf::from_file(paf);
            println!("{}", paf.sam_header());
            for rec in paf.records {
                println!("{}", rec.to_sam_string(fasta_reader.as_ref()));
            }
        }
        //
        // Run Fastx-split
        //
        Some(Commands::FastxSplit { fastx }) => {
            //run_split_fasta(fasta);
            fastx::run_split_fastx(fastx, "-");
        }
        //
        // Run GetFasta
        //
        Some(Commands::GetFasta {
            fasta,
            bed,
            strand,
            name,
        }) => {
            getfasta::get_fasta(fasta, bed, *name, *strand);
        }
        //
        // no command opt
        //
        None => {}
    };

    let duration = pg_start.elapsed();
    log::info!(
        "{} done! Time elapsed: {}",
        subcommand.bright_green().bold(),
        format!("{:.2?}", duration).bright_yellow().bold()
    );
}
