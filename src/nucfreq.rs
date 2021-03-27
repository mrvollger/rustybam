use regex::Regex;
use rust_htslib::{bam, bam::Read};
use std::fmt;
use std::str;

/// Code to get the coverage
/// # Example
/// ```
/// let rtn = rustybam::nucfreq::coverage("test/test.bam");
/// assert_eq!(rtn, 197);
/// ```
pub fn coverage(path: &str) -> u64 {
    eprintln!("Reading from {}", path);
    let mut bam = bam::Reader::from_path(path).unwrap();
    let mut count: u64 = 0;

    for p in bam.pileup() {
        let pileup = p.unwrap();
        //println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                let bp = alignment.record().seq()[alignment.qpos().unwrap()];
                if bp == b'A' {
                    count += 1;
                }
            }
        }
    }
    count
}

pub struct Nucfreq {
    pub pos: u32,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
}

impl fmt::Display for Nucfreq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.pos,
            self.pos + 1,
            self.a,
            self.c,
            self.g,
            self.t
        )
    }
}

/// get count for A,C,G,T at every pos
/// # Example
/// ```
/// use rust_htslib::{bam, bam::Read};
/// let mut bam = bam::IndexedReader::from_path("test/test_nucfreq.bam").unwrap();
/// let vec = rustybam::nucfreq::nucfreq( &mut bam );
/// eprintln!("{:?}", vec[0].a);
/// for f in vec {
///   let t = vec![f.a, f.c, f.g, f.t];
///   let max = t.iter().max().unwrap();
///   if(*max != 0){
///     assert_eq!(*max, 2);
///   }
/// }
/// ```
pub fn nucfreq(bam: &mut rust_htslib::bam::IndexedReader) -> Vec<Nucfreq> {
    let mut vec = Vec::new();
    for (_idx, p) in bam.pileup().enumerate() {
        let pileup = p.unwrap();
        //println!("{}", _idx);
        let mut freqs = Nucfreq {
            pos: pileup.pos(),
            a: 0,
            c: 0,
            g: 0,
            t: 0,
        };
        for aln in pileup.alignments() {
            if !aln.is_del() && !aln.is_refskip() {
                let bp = aln.record().seq()[aln.qpos().unwrap()];
                match bp {
                    b'A' => freqs.a += 1,
                    b'C' => freqs.c += 1,
                    b'G' => freqs.g += 1,
                    b'T' => freqs.t += 1,
                    b'N' => (),
                    _ => eprintln!("Seq character not recognized (u8):{:?}", bp),
                }
            }
        }
        vec.push(freqs);
    }
    vec
}

pub struct Region {
    pub name: String,
    pub st: u32,
    pub en: u32,
}

/// parse region strings
/// # Example
/// ```
/// let rgn = rustybam::nucfreq::parse_region("chr1:1-1000");
/// assert_eq!("chr1", rgn.name);
/// assert_eq!(0, rgn.st);
/// assert_eq!(1000, rgn.en);
///
/// let rgn2 = rustybam::nucfreq::parse_region("chr1:2-2000:1-1000");
/// assert_eq!("chr1:2-2000", rgn2.name);
/// ```
pub fn parse_region(region: &str) -> Region {
    let re = Regex::new(r"(.+):([0-9]+)-([0-9]+)").unwrap();
    let caps = re.captures(region).expect("Failed to parse region string.");

    Region {
        name: caps.get(1).unwrap().as_str().to_string(),
        st: caps.get(2).unwrap().as_str().parse::<u32>().unwrap() - 1,
        en: caps.get(3).unwrap().as_str().parse().unwrap_or(4294967295), // this is 2^32-1
    }
}

/// get count for A,C,G,T at every pos in the region
/// # Example
/// ```
/// use rust_htslib::{bam, bam::Read};
/// let mut bam = bam::IndexedReader::from_path("test/asm_small.bam").unwrap();
/// let vec = rustybam::nucfreq::region_nucfreq( &mut bam, "chr1:1-1000");
/// let vec2 = rustybam::nucfreq::region_nucfreq( &mut bam, "chr6:8-8000");
/// let vec3 = rustybam::nucfreq::region_nucfreq( &mut bam, "chr1:2-2000");
/// ```
pub fn region_nucfreq(bam: &mut rust_htslib::bam::IndexedReader, region: &str) -> Vec<Nucfreq> {
    let rgn = parse_region(region);
    eprintln!("Finding nucfreq in: {}\t{}\t{}", rgn.name, rgn.st, rgn.en);
    bam.fetch((&rgn.name, rgn.st as i64, rgn.en as i64))
        .unwrap();

    // get the nucfreq and filter for valid regions
    nucfreq(bam)
        .into_iter()
        .filter(|nf| nf.pos >= rgn.st && nf.pos < rgn.en)
        .collect()
}

pub fn print_nucfreq(vec: Vec<Nucfreq>) {
    for nf in vec {
        print!("{}", nf);
    }
}
