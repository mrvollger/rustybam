use bio::io::fasta;
use bio::io::fastq;
use itertools::Itertools;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustybam::cli::Commands;
use rustybam::paf::PafRecord;
use rustybam::*;
use std::collections::{HashMap, HashSet};
use std::io;
use std::time::Instant;

fn main() {
    parse_cli();
}

pub fn parse_cli() {
    let pg_start = Instant::now();
    let args = cli::make_cli_parse();

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
        Some(Commands::Bedlength { bed, readable }) => {
            let rgns = bed::parse_bed(bed);
            let count: u64 = rgns.into_iter().map(|rgn| rgn.en - rgn.st).sum();
            if *readable {
                println!("{}", (count as f64) / 1e6);
            } else {
                println!("{}", count);
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
        // Run Liftover
        //
        Some(Commands::Filter {
            paf,
            paired_len,
            aln,
            query,
        }) => {
            let mut paf = paf::Paf::from_file(paf);
            eprintln!("{} PAF records BEFORE filtering.", paf.records.len());
            paf.filter_query_len(*query);
            paf.filter_aln_len(*aln);
            paf.filter_aln_pairs(*paired_len);
            eprintln!("{} PAF records AFTER filtering.", paf.records.len());
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
            orient_records(paf, *scaffold, *insert);
        }
        //
        // Run Breakpaf
        //
        Some(Commands::Breakpaf { paf, max_size }) => {
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
        // Run Fasta-split
        //
        Some(Commands::FastaSplit { fasta }) => {
            run_split_fasta(fasta);
        }
        //
        // Run Fastq-split
        //
        Some(Commands::FastqSplit { fastq }) => {
            run_split_fastq(fastq);
        }
        //
        // Run Getfasta
        //
        Some(Commands::Getfasta {
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
    eprintln!("Time elapsed during rustybam: {:.3?}", duration);
}

pub fn run_split_fasta(files: &[String]) {
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
}

pub fn run_split_fastq(files: &[String]) {
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
}

pub fn orient_records(paf: &str, scaffold: bool, insert: u64) {
    let mut paf = paf::Paf::from_file(paf);
    let mut orient_order_dict = HashMap::new();
    let mut t_names = HashSet::new();
    // calculate whether a contig is mostly forward or reverse strand
    // and determine the middle alignment position with respect to the target
    for rec in &paf.records {
        let (orient, total_bp, order) = orient_order_dict
            .entry((rec.t_name.clone(), rec.q_name.clone()))
            .or_insert((0_i64, 0_u64, 0_u64));
        // set the orientation of the query relative to the target
        if rec.strand == '-' {
            *orient -= (rec.q_en - rec.q_st) as i64;
        } else {
            *orient += (rec.q_en - rec.q_st) as i64;
        }
        // set a number that will determine the order of the contig
        let weight = rec.t_en - rec.t_st;
        *total_bp += weight;
        *order += weight * (rec.t_st + rec.t_en) / 2;
        // make a list of targets
        t_names.insert(rec.t_name.clone());
    }

    // set the order and orientation of records
    for rec in &mut paf.records {
        // set the order of the records
        let (orient, total_bp, order) = orient_order_dict
            .get(&(rec.t_name.clone(), rec.q_name.clone()))
            .unwrap();
        rec.order = *order / *total_bp;

        // reverse record if it is mostly on the rc
        if *orient < 0 {
            rec.q_name = format!("{}-", rec.q_name);
            let new_st = rec.q_len - rec.q_en;
            let new_en = rec.q_len - rec.q_st;
            rec.q_st = new_st;
            rec.q_en = new_en;
            rec.strand = if rec.strand == '+' { '-' } else { '+' };
        } else {
            rec.q_name = format!("{}+", rec.q_name);
        }
    }

    // sort the records by their target name and order
    paf.records.sort_by(|a, b| {
        a.t_name
            .cmp(&b.t_name) // group by target
            .then(a.order.cmp(&b.order)) // order query by position in target
            .then(a.q_st.cmp(&b.q_st)) // order by position in query
    });

    // group by t_name
    for (_t_name, t_recs) in &paf.records.iter_mut().group_by(|rec| rec.t_name.clone()) {
        let mut t_recs: Vec<&mut PafRecord> = t_recs.collect();
        // sort recs by order
        t_recs.sort_by(|a, b| {
            a.order
                .cmp(&b.order) // order query by position in target
                .then(a.q_st.cmp(&b.q_st)) // order by position in query
        });

        // new scaffold name
        let scaffold_name = t_recs
            .iter()
            .map(|rec| rec.q_name.clone())
            .unique()
            .collect::<Vec<String>>()
            .join("::");

        let mut scaffold_len = 0_u64;
        for (_q_name, q_recs) in &t_recs.iter_mut().group_by(|rec| rec.q_name.clone()) {
            let q_recs: Vec<&mut &mut PafRecord> = q_recs.collect();
            let q_min = q_recs.iter().map(|rec| rec.q_st).min().unwrap_or(0);
            let q_max = q_recs.iter().map(|rec| rec.q_en).max().unwrap_or(0);
            let added_q_bases = q_max - q_min;
            for rec in q_recs {
                if scaffold {
                    rec.q_st = rec.q_st - q_min + scaffold_len;
                    rec.q_en = rec.q_en - q_min + scaffold_len;
                }
            }
            scaffold_len += added_q_bases + insert;
        }
        // remove padding insert on the end of rec
        scaffold_len -= insert;

        for rec in t_recs {
            if scaffold {
                rec.q_name = scaffold_name.clone();
                rec.q_len = scaffold_len;
            }
            // return result
            println!("{}", rec);
        }
    }
}
