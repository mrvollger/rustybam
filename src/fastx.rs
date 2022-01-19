use super::myio;
use needletail::{parse_fastx_file, parse_fastx_stdin, parser::LineEnding};

/// Split a fasta file across outputs
/// ```
/// use rustybam::fastx;
/// fastx::run_split_fastx(&["-".to_string()], ".test/large.test.fa.gz")
/// ```
pub fn run_split_fastx(files: &[String], infile: &str) {
    // open the output files
    let mut outs = Vec::new();
    for f in files {
        let handle = myio::writer(f);
        outs.push(handle);
    }
    // open reader
    let mut reader = if infile == "-" {
        parse_fastx_stdin().expect("Missing or invalid stdin for fastx parser.")
    } else {
        parse_fastx_file(infile).expect("Missing or invalid stdin for fastx parser.")
    };
    // iterate
    let mut out_idx = 0;
    let mut rec_num = 0;
    while let Some(record) = reader.next() {
        let seq_rec =
            record.unwrap_or_else(|_| panic!("Error reading record number {}", rec_num + 1));
        seq_rec
            .write(&mut outs[out_idx], Some(LineEnding::Unix))
            .unwrap_or_else(|_| panic!("Error writing record number {}", rec_num + 1));

        log::debug!("Wrote record number {}", rec_num + 1);
        out_idx += 1;
        rec_num += 1;
        if out_idx == outs.len() {
            out_idx = 0;
        }
    }
    // Close all the files.
    // this breaks things?
    for mut out in outs {
        out.flush()
            .unwrap_or_else(|_| panic!("Error flushing output!"));
    }
}
