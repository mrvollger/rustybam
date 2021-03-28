use super::bed::*;
use rust_htslib::bam::Read;
use std::fmt;

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

/// get count for A,C,G,T at every pos in the region
/// # Example
/// ```
/// use rust_htslib::{bam, bam::Read};
/// use rustybam::nucfreq::*;
/// use rustybam::bed::*;
///
/// let mut bam = bam::IndexedReader::from_path("test/asm_small.bam").unwrap();
///
/// let vec  = region_nucfreq( &mut bam, &parse_region("chr22:1-1000"));
/// let vec2 = region_nucfreq( &mut bam, &parse_region("chr21:8-8000"));
/// let vec3 = region_nucfreq( &mut bam, &parse_region("chr20:2-2000"));
/// ```
pub fn region_nucfreq(bam: &mut rust_htslib::bam::IndexedReader, rgn: &Region) -> Vec<Nucfreq> {
    eprintln!("Finding nucfreq in: {}\t{}\t{}", rgn.name, rgn.st, rgn.en);
    bam.fetch((&rgn.name, rgn.st as i64, rgn.en as i64))
        .unwrap();

    // get the nucfreq and filter for valid regions
    nucfreq(bam)
        .into_iter()
        .filter(|nf| nf.pos >= rgn.st && nf.pos < rgn.en)
        .collect()
}

pub fn print_nucfreq_header() {
    print!("#chr\tstart\tend\t");
    print!("A\tC\tG\tT\t");
    println!("region_id");
}

pub fn print_nucfreq(vec: Vec<Nucfreq>, rgn: &Region) {
    for nf in vec {
        print!("{}\t{}\t{}\t", rgn.name, nf.pos, nf.pos + 1);
        print!("{}\t{}\t{}\t{}\t", nf.a, nf.c, nf.g, nf.t);
        println!("{}", rgn.id);
    }
}
