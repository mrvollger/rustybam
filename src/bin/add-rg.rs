use clap::{App, Arg};
use log;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{self, Header, Read};

fn main() {
    let matches = App::new("add_rg")
        .version("0.1")
        .author("Your Name <your.email@example.com>")
        .about("Adds RG lines from a provided BAM file to a BAM file from stdin to stdout")
        .arg(
            Arg::new("source")
                .help("Source BAM file to read RG lines from")
                .required(true)
                .index(1),
        )
        // add threads argument with default of 8
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use")
                .takes_value(true)
                .default_value("8"),
        )
        .get_matches();

    // set log level to info
    env_logger::Builder::new()
        .filter(None, log::LevelFilter::Info)
        .init();

    let n_threads = matches
        .value_of("threads")
        .unwrap_or("8")
        .parse::<usize>()
        .unwrap();
    let source_file = matches.value_of("source").unwrap();

    // Open the source BAM file and read its header
    let source_bam = bam::Reader::from_path(source_file).expect("Failed to open source BAM file");
    let source_header = source_bam.header();
    let source_header_text = String::from_utf8_lossy(source_header.as_bytes());

    // Extract RG lines from the source header
    let rg_lines: Vec<&str> = source_header_text
        .lines()
        .filter(|line| line.starts_with("@RG"))
        // strip the leading @ character
        .map(|line| line.trim_start_matches('@'))
        .collect();

    if rg_lines.is_empty() {
        log::warn!("No RG lines found in the source BAM file.");
    }

    // Open the target BAM file and read its header
    let mut target_bam = bam::Reader::from_stdin().expect("Failed to open BAM file from stdin");
    target_bam
        .set_threads(n_threads)
        .expect("Failed to set threads");

    let target_header = target_bam.header();
    let mut new_header = Header::from_template(target_header);

    // Add RG lines to the new header
    for rg_line in rg_lines {
        let header_line = HeaderRecord::new(rg_line.as_bytes());
        new_header.push_record(&header_line);
    }

    // Create a writer for the output BAM file
    let mut output_bam = bam::Writer::from_stdout(&new_header, bam::Format::Bam)
        .expect("Failed to create output BAM writer");
    output_bam
        .set_threads(n_threads)
        .expect("Failed to set threads for output BAM writer");

    // Write records from the target BAM file to the output BAM file
    for record in target_bam.records() {
        let record = record.expect("Failed to read record from target BAM file");
        output_bam
            .write(&record)
            .expect("Failed to write record to output BAM file");
    }

    println!("RG lines successfully added to the output BAM file.");
}
