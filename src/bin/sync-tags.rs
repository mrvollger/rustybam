use clap::Parser;
use log;
use rust_htslib::bam::{self, Read};
use std::fmt::Debug;

// Args,
#[derive(Parser, Debug)]
/// This program synchronizes tags between two BAM files. It applies any tags from the first BAM file to the second BAM file if the read names match and the tags do not already exist in the second BAM file and writes the output to stdout. The order of the reads must be the same in both BAM files. This can be done with, e.g. `samtools sort -N`.
pub struct SyncTagsArgs {
    /// First BAM file (source of tags)
    bam1: String,
    /// Second BAM file (tags will be updated)
    bam2: String,
    /// Output BAM file
    #[clap(short, long, default_value = "-")]
    output: String,
    /// Number of threads to use
    #[clap(short, long, default_value_t = 8)]
    threads: usize,
    /// Write uncompressed output
    #[clap(short = 'u', long)]
    uncompressed: bool,
}

pub fn bam_reader_from_path_or_stdin(path: &str, threads: usize) -> bam::Reader {
    let mut bam = if path == "-" {
        bam::Reader::from_stdin().expect("Failed to open BAM file from stdin")
    } else {
        bam::Reader::from_path(path).expect("Failed to open BAM file from path")
    };
    bam.set_threads(threads)
        .expect("Failed to set threads for BAM reader");
    bam
}

pub fn bam_writer_from_path_or_stdout(
    path: &str,
    header: &bam::HeaderView,
    threads: usize,
    uncompressed: bool,
) -> bam::Writer {
    let mut header = bam::Header::from_template(header);
    // add a PG line to the header
    let mut pg_line = bam::header::HeaderRecord::new(b"PG");
    pg_line.push_tag(b"ID", "sync-tags");
    pg_line.push_tag(b"PN", "sync-tags");
    pg_line.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
    // get the full command line call as a string
    let full_cmd = std::env::args()
        .map(|arg| arg.replace(' ', "\\ "))
        .collect::<Vec<String>>()
        .join(" ");
    pg_line.push_tag(b"CL", full_cmd);
    header.push_record(&pg_line);

    let mut writer = if path == "-" {
        bam::Writer::from_stdout(&header, bam::Format::Bam)
            .expect("Failed to create BAM writer for stdout")
    } else {
        bam::Writer::from_path(path, &header, bam::Format::Bam)
            .expect("Failed to create BAM writer from path")
    };
    writer
        .set_threads(threads)
        .expect("Failed to set threads for BAM writer");
    if uncompressed {
        writer
            .set_compression_level(bam::CompressionLevel::Uncompressed)
            .expect("Failed to set compression level");
    }
    writer
}

fn main() {
    let args = SyncTagsArgs::parse();

    // Set log level to info
    env_logger::Builder::new()
        .filter(None, log::LevelFilter::Info)
        .init();

    // Open the first BAM file (source of tags)
    let mut bam1 = bam_reader_from_path_or_stdin(&args.bam1, args.threads);

    // Open the second BAM file (tags will be updated)
    let mut bam2 = bam_reader_from_path_or_stdin(&args.bam2, args.threads);

    // Create a writer for the output BAM file
    let header = bam2.header();
    let mut output_bam =
        bam_writer_from_path_or_stdout(&args.output, header, args.threads, args.uncompressed);

    let bam1_iter = bam1.records();
    let mut bam2_iter = bam2.records();
    let mut destination_rec = if let Some(destination_rec) = bam2_iter.next() {
        destination_rec.expect("Failed to read record from second BAM file")
    } else {
        log::warn!("No records in the second BAM file.");
        return;
    };

    // Iterate over the first BAM file
    for template_rec in bam1_iter {
        let template_rec = template_rec.expect("Failed to read record from first BAM file");
        // Iterate over the second BAM file
        while template_rec.qname() == destination_rec.qname() {
            // Move the tags from the template record to the next read
            for (key, value) in template_rec
                .aux_iter()
                .map(|x| x.expect("Unable to read tags from template BAM file"))
            {
                if !destination_rec.aux(key).is_ok() {
                    destination_rec
                        .push_aux(key, value)
                        .expect("Failed to push tag to next read");
                }
            }

            // Write the updated read to the output BAM file
            output_bam
                .write(&destination_rec)
                .expect("Failed to write record to output BAM file");

            // get the next read
            destination_rec = if let Some(next_read) = bam2_iter.next() {
                next_read.expect("Failed to read record from second BAM file")
            } else {
                log::warn!("No more records in the second BAM file.");
                break;
            }
        }
    }

    log::info!("Tags successfully synchronized and written to output BAM file.");
}
