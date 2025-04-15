use log;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{self, Header, Read};

pub fn add_rg(threads: usize, source_file: &str, uncompressed: bool) {
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
        // keep only unique lines
        .collect::<std::collections::HashSet<_>>()
        // convert to Vec
        .into_iter()
        .collect();

    if rg_lines.is_empty() {
        log::warn!("No RG lines found in the source BAM file.");
    }

    // Open the target BAM file and read its header
    let mut target_bam = bam::Reader::from_stdin().expect("Failed to open BAM file from stdin");
    target_bam
        .set_threads(threads)
        .expect("Failed to set threads");

    let target_header = target_bam.header();
    let mut new_header = Header::from_template(target_header);

    // Add RG lines to the new header
    for rg_line in rg_lines {
        let header_line = HeaderRecord::new(rg_line.as_bytes());
        new_header.push_record(&header_line);
    }
    // remove duplicate RG lines

    // Create a writer for the output BAM file
    let mut output_bam = bam::Writer::from_stdout(&new_header, bam::Format::Bam)
        .expect("Failed to create output BAM writer");
    output_bam
        .set_threads(threads)
        .expect("Failed to set threads for output BAM writer");

    // Check if the output should be uncompressed
    if uncompressed {
        output_bam
            .set_compression_level(bam::CompressionLevel::Uncompressed)
            .expect("Failed to set compression level for output BAM writer");
    }

    // Write records from the target BAM file to the output BAM file
    for record in target_bam.records() {
        let record = record.expect("Failed to read record from target BAM file");

        output_bam
            .write(&record)
            .expect("Failed to write record to output BAM file");
    }
    log::info!("RG lines successfully added to the output BAM file.");
}
