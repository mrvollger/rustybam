use super::bed::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::convert::TryFrom;
use std::fmt;

pub struct Nucfreq {
    pub name: String,
    pub pos: u32,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub id: String,
}

impl fmt::Display for Nucfreq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.name,
            self.pos,
            self.pos + 1,
            self.a,
            self.c,
            self.g,
            self.t,
            self.id
        )
    }
}

pub struct NucList {
    pub rgn: Region,
    pub nucfreqs: Vec<Nucfreq>,
}

/// get count for A,C,G,T at every pos
/// # Example
/// ```
/// use rust_htslib::{bam, bam::Read};
/// let mut bam = bam::IndexedReader::from_path(".test/test_nucfreq.bam").unwrap();
/// let rgn = rustybam::bed::Region {
///     name : "CHROMOSOME_I".to_string(),
///     st :  1,
///     en : 102,
///     id : "None".to_string(),
///     ..Default::default()
/// };
/// let vec = rustybam::nucfreq::nucfreq( &mut bam, &rgn);
/// eprintln!("{:?}", vec[0].a);
/// for f in vec {
///   let t = vec![f.a, f.c, f.g, f.t];
///   let max = t.iter().max().unwrap();
///   if(*max != 0){
///     assert_eq!(*max, 2);
///   }
/// }
/// ```
pub fn nucfreq(bam: &mut rust_htslib::bam::IndexedReader, rgn: &Region) -> Vec<Nucfreq> {
    let rgn_length = usize::try_from(rgn.en - rgn.st + 1).unwrap();
    let mut vec = Vec::with_capacity(rgn_length);
    for p in bam.pileup() {
        let pileup = p.unwrap();
        if (pileup.pos() as u64) < rgn.st || (pileup.pos() as u64) >= rgn.en {
            continue;
        }
        //println!("{}", _idx);
        let mut freqs = Nucfreq {
            name: rgn.name.clone(),
            pos: pileup.pos(),
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            id: rgn.id.clone(),
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
/// let mut bam = ".test/asm_small.bam";
///
/// let vec  = region_nucfreq( bam, &parse_region("chr22:1-1000"), 1);
/// let vec2 = region_nucfreq( bam, &parse_region("chr21:8-8000"), 4);
/// let vec3 = region_nucfreq( bam, &parse_region("chr20:2-2000"), 2);
/// ```
pub fn region_nucfreq(bam_f: &str, rgn: &Region, threads: usize) -> Vec<Nucfreq> {
    eprint!("\rFinding nucfreq in: {}\t{}\t{}", rgn.name, rgn.st, rgn.en);
    // open the bam
    let mut bam =
        bam::IndexedReader::from_path(bam_f).unwrap_or_else(|_| panic!("Failed to open {}", bam_f));
    bam.set_threads(threads)
        .expect("Failed to enable multithreaded bam reading.");

    // read the bam
    bam.fetch((&rgn.name, rgn.st as i64, rgn.en as i64))
        .unwrap_or_else(|_| panic!("Is this region ({}) in your reference/bam?", rgn));

    // get the nucfreq and filter for valid regions
    nucfreq(&mut bam, rgn)
}

pub fn print_nucfreq_header() {
    print!("#chr\tstart\tend\t");
    print!("A\tC\tG\tT\t");
    println!("region_id");
}

pub fn print_nucfreq(vec: &[Nucfreq]) {
    for nf in vec {
        println!("{}", nf);
    }
}

pub fn small_nucfreq(vec: &[Nucfreq]) {
    let mut cur_name = "".to_string();
    let mut cur_id = "".to_string();

    for nf in vec {
        if nf.name != cur_name || nf.id != cur_id {
            cur_name = nf.name.clone();
            cur_id = nf.id.clone();
            println!("#{}\t{}\t{}", nf.name, nf.pos, nf.id);
        }
        let mut mc = [nf.a, nf.c, nf.g, nf.t];
        mc.sort_unstable();
        println!("{}\t{}", mc[3], mc[2]);
    }
}
