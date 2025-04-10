use clap::{App, Arg};
use niffler::get_reader;
use num_format::{Locale, ToFormattedString};
use rust_htslib::bam::{self, Read};
use std::fs;
use std::io::{self, BufRead};

fn read_bam(file: &str, threads: usize) -> Option<(Vec<String>, Vec<usize>)> {
    let mut names = Vec::new();
    let mut lengths = Vec::new();

    let mut bam = bam::Reader::from_path(file).ok()?;
    bam.set_threads(threads).ok()?;
    for record in bam.records() {
        let rec = record.ok()?;
        if rec.is_unmapped() || (!rec.is_secondary() && !rec.is_supplementary()) {
            names.push(String::from_utf8_lossy(rec.qname()).to_string());
            lengths.push(rec.seq().len());
        }
    }
    eprintln!("SAM/BAM read: {}", file);
    Some((names, lengths))
}

fn read_bed(file: &str) -> Option<(Vec<String>, Vec<usize>)> {
    let mut names = Vec::new();
    let mut lengths = Vec::new();

    let file = fs::File::open(file).ok()?;
    let (reader, _compression) = get_reader(Box::new(file)).ok()?;
    let reader = io::BufReader::new(reader);

    for line in reader.lines() {
        let line = line.ok()?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 3 {
            let start: usize = fields[1].parse().ok()?;
            let end: usize = fields[2].parse().ok()?;
            names.push(fields[0].to_string());
            lengths.push(end - start);
        }
    }
    Some((names, lengths))
}

fn calc_stats(
    lengths: &[usize],
    quantiles: &[f64],
    genome_size: Option<usize>,
) -> (usize, usize, f64, Vec<f64>, usize, usize, usize, f64) {
    let n = lengths.len();
    let total: usize = genome_size.unwrap_or_else(|| lengths.iter().sum());
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));

    let max = *sorted_lengths.first().unwrap_or(&0);
    let min = *sorted_lengths.last().unwrap_or(&0);
    let mean = total as f64 / n as f64;

    let au_n: f64 = sorted_lengths.iter().map(|&x| (x * x) as f64).sum::<f64>() / total as f64;

    let mut quantile_values = Vec::new();
    for &q in quantiles {
        let idx = (q * n as f64).ceil() as usize - 1;
        quantile_values.push(sorted_lengths.get(idx).cloned().unwrap_or(0) as f64);
    }

    let mut cumulative = 0;
    let mut n50 = 0;
    for &len in &sorted_lengths {
        cumulative += len;
        if cumulative >= total / 2 {
            n50 = len;
            break;
        }
    }

    (total, n, mean, quantile_values, min, max, n50, au_n)
}

pub fn h_fmt<T>(num: T) -> String
where
    T: Into<f64> + Copy,
{
    let mut num: f64 = num.into();
    for unit in ["", "Kbp", "Mbp"] {
        if num < 1000.0 {
            return format!("{:.2}{}", num, unit);
        }
        num /= 1000.0;
    }
    format!("{:.2}{}", num, "Gbp")
}

fn main() {
    let matches = App::new("seq_stats")
        .version("0.1")
        .author("Mitchell R. Vollger <mrvollger@gmail.com>")
        .about("Calculates sequence statistics from various file formats")
        .arg(
            Arg::new("infiles")
                .help("Input files (fast{a,q}(.gz), sam, bam, bed)")
                .required(true)
                .multiple_values(true),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use")
                .takes_value(true)
                .default_value("4"),
        )
        .arg(
            Arg::new("human")
                .short('r')
                .long("human")
                .help("Print human-readable output")
                .takes_value(false),
        )
        .arg(
            Arg::new("quantiles")
                .short('q')
                .long("quantiles")
                .help("Quantiles to calculate")
                .takes_value(true)
                .multiple_values(true)
                .default_value("0.5"),
        )
        .arg(
            Arg::new("genome_size")
                .short('g')
                .long("genome_size")
                .help("Genome size for NG50 calculation")
                .takes_value(true),
        )
        .get_matches();

    // set log level to info
    env_logger::Builder::new()
        .filter(None, log::LevelFilter::Info)
        .init();

    let infiles: Vec<&str> = matches.values_of("infiles").unwrap().collect();
    let threads: usize = matches.value_of_t("threads").unwrap_or(4);
    let human_readable = matches.is_present("human");
    let quantiles: Vec<f64> = matches
        .values_of("quantiles")
        .unwrap()
        .map(|q| q.parse().unwrap_or(0.5))
        .collect();
    let genome_size: Option<usize> = matches.value_of_t("genome_size").ok();

    let mut output = String::from("file\ttotalBp\tnSeqs\tmean\tquantiles\tmin\tmax\tN50\tauN\n");

    for file in infiles {
        let lengths = if file.ends_with(".bam") || file.ends_with(".sam") || file.ends_with(".cram")
        {
            log::info!("Reading BAM/SAM/CRAM file: {}", file);
            read_bam(file, threads)
        } else if file.ends_with(".bed") || file.ends_with(".bed.gz") {
            log::info!("Reading BED file: {}", file);
            read_bed(file)
        } else {
            None
        };

        if let Some((_, lengths)) = lengths {
            let (total, n, mean, quantile_values, min, max, n50, au_n) =
                calc_stats(&lengths, &quantiles, genome_size);

            let quantile_str = quantile_values
                .iter()
                .map(|q| {
                    if human_readable {
                        h_fmt(*q)
                    } else {
                        q.to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join("\t");

            let line = if human_readable {
                format!(
                    "{}\t{}\t{}\t{:}\t{}\t{}\t{}\t{}\t{:}\n",
                    file,
                    h_fmt(total as f64),
                    n.to_formatted_string(&Locale::en),
                    h_fmt(mean),
                    quantile_str,
                    h_fmt(min as f64),
                    h_fmt(max as f64),
                    h_fmt(n50 as f64),
                    h_fmt(au_n)
                )
            } else {
                format!(
                    "{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\n",
                    file, total, n, mean, quantile_str, min, max, n50, au_n
                )
            };
            output.push_str(&line);
        } else {
            eprintln!("Skipping file: {}", file);
        }
    }

    print!("{}", output);
}
