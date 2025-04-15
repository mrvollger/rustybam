use linear_map::LinearMap;
use log;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{self, Header, Read};
use std::collections::HashMap;

pub fn header_from_hashmap(hash_header: HashMap<String, Vec<LinearMap<String, String>>>) -> Header {
    let mut header = Header::new();
    for (key, values) in hash_header.iter() {
        for value in values {
            let mut record = HeaderRecord::new(key.as_bytes());
            for (tag, val) in value.iter() {
                record.push_tag(tag.as_bytes(), val);
            }
            header.push_record(&record);
        }
    }
    header
}

pub fn get_rg_ids(header: &Header) -> Vec<String> {
    let mut rg_ids = Vec::new();
    let hash = header.to_hashmap();
    for key in hash.keys() {
        if key.eq("RG") {
            for rg_line in hash.get(key).unwrap() {
                if let Some(id) = rg_line.get("ID") {
                    rg_ids.push(id.to_string());
                }
            }
        }
    }
    rg_ids
}

pub fn add_rg(threads: usize, source_file: &str, uncompressed: bool, sample: &Option<String>) {
    // Open the source BAM file and read its header
    let source_bam = bam::Reader::from_path(source_file).expect("Failed to open source BAM file");
    let source_header_view = source_bam.header();
    let source_header = Header::from_template(source_header_view);
    let mut source_header_hash = source_header.to_hashmap();

    // Extract RG lines from the source header
    let rg_ids_added = get_rg_ids(&source_header);
    if rg_ids_added.is_empty() {
        log::warn!("No RG lines found in the source BAM file. None will be added.");
    }

    // Open the target BAM file and read its header
    let mut target_bam = bam::Reader::from_stdin().expect("Failed to open BAM file from stdin");
    target_bam
        .set_threads(threads)
        .expect("Failed to set threads");

    let target_header = target_bam.header();
    let mut new_header = Header::from_template(target_header).to_hashmap();

    // Add RG lines from the source header to the new header
    let new_header_rg_lines = new_header.entry("RG".to_string()).or_default();

    // remove rg lines with IDs in the target header that are also in the source header
    let fake_id = "".to_string();
    new_header_rg_lines.retain(|rg_line| {
        let id = rg_line.get("ID").unwrap_or(&fake_id);
        !rg_ids_added.contains(id)
    });

    // get RG lines from the source header
    let source_rg_lines = source_header_hash.entry("RG".to_string()).or_default();

    // Add sample if provided to the source rg lines
    if let Some(sample) = sample {
        for rg_line in source_rg_lines.iter_mut() {
            rg_line.insert("SM".to_string(), sample.clone());
        }
    }

    // Add RG lines from the source header to the new header
    new_header_rg_lines.extend(source_rg_lines.iter().cloned());

    // Create a writer for the output BAM file
    let out_header = header_from_hashmap(new_header);
    let mut output_bam = bam::Writer::from_stdout(&out_header, bam::Format::Bam)
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
