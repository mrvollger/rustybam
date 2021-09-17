use flate2::read;
/*use flate2::read::MultiGzDecoder;
use flate2::write;
use flate2::write::DeflateDecoder;*/
use flate2::Compression;
use gzp::deflate::Bgzf; //, Gzip, Mgzip, RawDeflate};
use gzp::ZBuilder;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const BUFFERSIZE: usize = 64 * 1024;

/// code modified from https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561/2
/// to work with gzp

/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &str) -> Box<dyn Write> {
    let path = Path::new(filename);
    let file = match File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        let buffer = Box::new(BufWriter::with_capacity(BUFFERSIZE, file));
        let writer = ZBuilder::<Bgzf, _>::new()
            .num_threads(8)
            .compression_level(Compression::new(9))
            .from_writer(buffer);
        Box::new(writer)
        //let writer = write::GzEncoder::new(file, Compression::default());
        //Box::new(BufWriter::with_capacity(128 * 1024, writer))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    }
}

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}
