use super::bed;
use core::fmt;
use regex::Regex;
use rust_htslib::bam::record::Cigar::*;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::record::*;
use std::io::BufRead;
//use std::convert::TryFrom;
use std::str::FromStr;

#[derive(Debug, Default, Clone)]
pub struct PafRecord {
    pub q_name: String,
    pub q_len: u64,
    pub q_st: u64,
    pub q_en: u64,
    pub strand: char,
    pub t_name: String,
    pub t_len: u64,
    pub t_st: u64,
    pub t_en: u64,
    pub nmatch: u64,
    pub aln_len: u64,
    pub mapq: u64,
    pub cigar: String,
    pub cs: String,
}

impl fmt::Display for PafRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.q_name,
            self.q_len,
            self.q_st,
            self.q_en,
            self.strand,
            self.t_name,
            self.t_len,
            self.t_st,
            self.t_en,
            self.nmatch,
            self.aln_len,
            self.mapq
        )
    }
}

pub fn consumes_reference(cigar_opt: &Cigar) -> bool {
    matches!(
        cigar_opt,
        Match(_i) | Del(_i) | RefSkip(_i) | Diff(_i) | Equal(_i)
    )
}
/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::Cigar::*;
/// assert!(paf::consumes_query(&Diff(5)));
/// ```
pub fn consumes_query(cigar_opt: &Cigar) -> bool {
    matches!(
        cigar_opt,
        Match(_i) | Ins(_i) | SoftClip(_i) | Diff(_i) | Equal(_i)
    )
}

pub fn paf_overlaps_target(paf: &PafRecord, rgn: &bed::Region) -> bool {
    if paf.t_name != rgn.name {
        return false;
    }
    paf.t_en > rgn.st && paf.t_st < rgn.en
}

/// Create a CigarString from given str.
/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::*;
/// use rust_htslib::bam::record::CigarString;
/// use rust_htslib::bam::record::Cigar::*;
/// use std::convert::TryFrom;
/// use std::str::FromStr;
/// let cigars = vec!["10M4D100I1102=", "100000M20=5P10X4M"];
/// for cigar_str in cigars{
///     let my_parse = paf::cigar_from_str(cigar_str).expect("Unable to parse cigar");
///     let hts_parse = CigarString::try_from(cigar_str).expect("Unable to parse cigar");
///     assert_eq!(my_parse, hts_parse);
/// }
/// ```
/// TODO add IO errors
pub fn cigar_from_str(text: &str) -> Result<CigarString, std::io::Error> {
    let bytes = text.as_bytes();
    let mut inner = Vec::new();
    let mut i = 0;
    let text_len = text.len();
    while i < text_len {
        let mut j = i;
        while j < text_len && bytes[j].is_ascii_digit() {
            j += 1;
        }
        let n = u32::from_str(&text[i..j]).expect("Bam parse error");
        let op = &text[j..j + 1];
        inner.push(match op {
            "M" => Cigar::Match(n),
            "I" => Cigar::Ins(n),
            "D" => Cigar::Del(n),
            "N" => Cigar::RefSkip(n),
            "H" => Cigar::HardClip(n),
            "S" => Cigar::SoftClip(n),
            "P" => Cigar::Pad(n),
            "=" => Cigar::Equal(n),
            "X" => Cigar::Diff(n),
            op => panic!("Unable to parse {} as cigar op", op),
        });
        i = j + 1;
    }
    Ok(CigarString(inner))
}

pub fn trim_paf_rec_to_rgn(rgn: &bed::Region, paf: &PafRecord) -> PafRecord {
    //let bytes = paf.cigar.as_bytes();
    //let cigar = CigarString::try_from(bytes).expect("Unable to parse cigar");
    let cigar = cigar_from_str(paf.cigar.as_str()).expect("Unable to parse cigar.");

    // initalize a trimmed paf record
    let mut trimmed_paf = (*paf).clone();

    // count until we get to the correct regions
    let mut t_pos = paf.t_st;
    let mut q_pos = paf.q_st;
    for opt in cigar.into_iter() {
        let moves_t = consumes_reference(&opt);
        let moves_q = consumes_query(&opt);
        let opt_len = opt.len() as u64;
        // incrment the ref and or query
        if moves_t {
            t_pos += opt_len;
        }
        if moves_q {
            q_pos += opt_len;
        }
        // record the starting positions
        if moves_t && t_pos - opt_len < rgn.st && t_pos >= rgn.st {
            trimmed_paf.t_st = rgn.st;
            trimmed_paf.q_st = q_pos;
            // correct q_st if we overshot with the last operation
            if moves_q {
                trimmed_paf.q_st -= t_pos - rgn.st;
            }
        }
        // record the ending positions
        if moves_t && t_pos >= rgn.en {
            trimmed_paf.t_en = rgn.en;
            trimmed_paf.q_en = q_pos;
            // correct q_en if we overshot with the last operation
            if moves_q {
                trimmed_paf.q_en -= t_pos - rgn.en;
            }
            break;
        }
    }
    // fix strand if needed
    if paf.strand == '-' {
        let q_en_tmp = trimmed_paf.q_en;
        trimmed_paf.q_en = trimmed_paf.q_len - trimmed_paf.q_st;
        trimmed_paf.q_st = trimmed_paf.q_len - q_en_tmp;
    }
    trimmed_paf
}

pub fn trim_paf_to_rgn(rgn: &bed::Region, paf: &[PafRecord]) -> Vec<PafRecord> {
    let mut trimmed_paf = Vec::new();
    //let mut i = 0;
    for rec in paf {
        //i += 1; eprintln!("{}", i);
        if paf_overlaps_target(rec, rgn) {
            trimmed_paf.push(trim_paf_rec_to_rgn(rgn, rec));
        }
    }
    trimmed_paf
}

/// # Example
/// ```
/// use rustybam::paf;
///
/// paf::read_paf_line("A 1 2 3 + B 1 2 3 10 11 60");
/// let rec = paf::make_fake_paf_rec();
/// assert_eq!("4M1I1D3=", rec.cigar);
///
/// ```
pub fn read_paf_line(line: &str) -> PafRecord {
    let t: Vec<&str> = line.split_ascii_whitespace().collect();
    assert!(t.len() >= 12); // must have all required columns
    let mut rec = PafRecord {
        q_name: t[0].to_string(),
        q_len: t[1].parse::<u64>().expect("q_len must be u64"),
        q_st: t[2].parse::<u64>().expect("q_st must be u64"),
        q_en: t[3].parse::<u64>().expect("q_en must be u64"),
        strand: t[4].parse::<char>().expect("strand must be a single char"),
        t_name: t[5].to_string(),
        t_len: t[6].parse::<u64>().expect("t_len must be u64"),
        t_st: t[7].parse::<u64>().expect("t_st must be u64"),
        t_en: t[8].parse::<u64>().expect("t_en must be u64"),
        nmatch: t[9].parse::<u64>().expect("nmatch must be u64"),
        aln_len: t[10].parse::<u64>().expect("aln_len must be u64"),
        mapq: t[11].parse::<u64>().expect("mapq must be u64"),
        ..Default::default()
    };

    let pattern = Regex::new("(..):(.):(.*)").unwrap();
    for token in t.iter().skip(12) {
        assert!(pattern.is_match(token));
        let caps = pattern.captures(token).unwrap();
        let tag = &caps[1];
        let value = &caps[3];
        if tag == "cg" {
            rec.cigar = value.to_string();
        } else if tag == "cs" {
            rec.cs = value.to_string();
        }
    }
    rec
}

/// # Example
/// ```
/// use rustybam::paf;
/// use std::fs::File;
/// use std::io::*;
/// let paf_file = Box::new(BufReader::new(File::open("test/asm_small.paf").expect("Unable to open file")));
/// let paf_recs = paf::read_paf(paf_file);
/// assert_eq!(paf_recs.len(), 249);
///
/// ```
pub fn read_paf(paf_file: Box<dyn BufRead>) -> Vec<PafRecord> {
    let mut vec = Vec::new();
    for (_index, line) in paf_file.lines().enumerate() {
        vec.push(read_paf_line(&line.unwrap()));
    }
    vec
}

pub fn make_fake_paf_rec() -> PafRecord {
    read_paf_line("Q 10 2 8 - T 20 12 18 3 8 60 cg:l:4M1I1D3=")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::Region;

    #[test]
    fn test_paf_trim() {
        let paf = make_fake_paf_rec();
        let rgn = Region {
            name: "T".to_string(),
            st: 14,
            en: 16,
            id: "None".to_string(),
        };
        let trim_paf = trim_paf_rec_to_rgn(&rgn, &paf);
        eprintln!("{:?}", trim_paf);
    }
}
