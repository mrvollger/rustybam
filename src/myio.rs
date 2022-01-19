use anyhow::Result;
use flate2::read;
use gzp::deflate::Bgzf; //, Gzip, Mgzip, RawDeflate};
use gzp::BgzfSyncReader;
use gzp::Compression;
use gzp::ZBuilder;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

type DynResult<T> = Result<T, Box<dyn Error + 'static>>;
const BUFFER_SIZE: usize = 32 * 1024;

/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &str) -> Box<dyn Write> {
    let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);
    let buffer = get_output(Some(path)).expect("Error: cannot create output file");

    if ext == Some(OsStr::new("gz")) {
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
/// Uses the presence of a `.gz` or `.bgz` extension to decide
/// Examples with zipped and unzipped
/// ```
/// use rustybam::myio::*;
/// use std::io::BufRead; // must import BufRead or you get an error at `lines()`
/// let n_paf = reader(".test/asm_small.paf").lines().count();
/// let n_paf_bgz = reader(".test/asm_small.paf.bgz").lines().count();
/// assert_eq!(n_paf, n_paf_bgz);
/// let n_paf_gz = reader(".test/asm_small.paf.gz").lines().count();
/// assert_eq!(n_paf, n_paf_gz);
/// ```
pub fn reader(filename: &str) -> Box<dyn BufRead + Send + 'static> {
    //Box<dyn BufRead> {
    let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);

    if ext == Some(OsStr::new("gz")) {
        let file = match File::open(&path) {
            Err(why) => panic!("couldn't open {}: {}", path.display(), why),
            Ok(file) => file,
        };
        Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            read::GzDecoder::new(file),
        ))
    } else if ext == Some(OsStr::new("bgz")) {
        Box::new(BufReader::new(BgzfSyncReader::new(
            get_input(Some(path)).expect("Error: cannot read input file."),
        )))
    } else {
        get_input(Some(path)).expect("Error: cannot read input file")
    }
}

/// Get a buffered output writer from stdout or a file
fn get_output(path: Option<PathBuf>) -> Result<Box<dyn Write + Send + 'static>> {
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

/// Get a buffered input reader from stdin or a file
fn get_input(path: Option<PathBuf>) -> DynResult<Box<dyn BufRead + Send + 'static>> {
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
