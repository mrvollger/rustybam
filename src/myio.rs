use anyhow::Result;
use flate2::read;
use flate2::Compression;
use gzp::deflate::Bgzf; //, Gzip, Mgzip, RawDeflate};
use gzp::ZBuilder;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::path::PathBuf;

const BUFFER_SIZE: usize = 128 * 1024;

/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &str) -> Box<dyn Write> {
    let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);
    let buffer = get_output(Some(path)).expect("Error: cannot create output file");

    if ext == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        let writer = ZBuilder::<Bgzf, _>::new()
            .num_threads(8)
            .compression_level(Compression::new(6))
            .from_writer(buffer);
        Box::new(writer)
    } else {
        buffer
    }
}

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);

    if ext == Some(OsStr::new("gz")) {
        let file = match File::open(&path) {
            Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
            Ok(file) => file,
        };
        Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            read::GzDecoder::new(file),
        ))
    } else {
        get_input(Some(path)).expect("Error: cannot read input file")
    }
}

/// Get a buffered output writer from stdout or a file
pub fn get_output(path: Option<PathBuf>) -> Result<Box<dyn Write + Send + 'static>> {
    let writer: Box<dyn Write + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, io::stdout()))
            } else {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, File::create(path)?))
            }
        }
        None => Box::new(BufWriter::with_capacity(BUFFER_SIZE, io::stdout())),
    };
    Ok(writer)
}

/// Get a bufferd input reader from stdin or a file
pub fn get_input(path: Option<PathBuf>) -> Result<Box<dyn BufRead + Send + 'static>> {
    let reader: Box<dyn BufRead + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, io::stdin()))
            } else {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, File::open(path)?))
            }
        }
        None => Box::new(BufReader::with_capacity(BUFFER_SIZE, io::stdin())),
    };
    Ok(reader)
}
