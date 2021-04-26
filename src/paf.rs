use super::bed;
use core::fmt;
use regex::Regex;
use rust_htslib::bam::record::Cigar::*;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::record::*;
use std::io::BufRead;
//use std::convert::TryFrom;
use std::str::FromStr;

#[derive(Debug)]
pub enum Error {
    PafParseCigar { msg: String },
    ParseIntError { msg: String },
    ParsePafColumn {},
}
type PafResult<T> = Result<T, crate::paf::Error>;

#[derive(Debug, Clone)]
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
    pub cigar: CigarString,
    pub tags: String,
}

impl fmt::Display for PafRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}{}",
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
            self.mapq,
            self.cigar.to_string(),
            self.tags
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
pub fn cigar_from_str(text: &str) -> PafResult<CigarString> {
    let bytes = text.as_bytes();
    let mut inner = Vec::new();
    let mut i = 0;
    let text_len = text.len();
    while i < text_len {
        let mut j = i;
        while j < text_len && bytes[j].is_ascii_digit() {
            j += 1;
        }
        let n = u32::from_str(&text[i..j]).map_err(|_| Error::PafParseCigar {
            msg: "expected integer".to_owned(),
        })?;
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
            op => {
                return Err(Error::PafParseCigar {
                    msg: format!("Cannot parse opt: {}", op),
                })
            }
        });
        i = j + 1;
    }
    Ok(CigarString(inner))
}

/// Swaps the query and reference and inverts the cigar sting
pub fn paf_swap_query_and_target(paf: &PafRecord) -> PafRecord {
    let mut flipped = paf.clone();
    // flip target
    flipped.t_name = paf.q_name.clone();
    flipped.t_len = paf.q_len;
    flipped.t_st = paf.q_st;
    flipped.t_en = paf.q_en;
    // flip query
    flipped.q_name = paf.t_name.clone();
    flipped.q_len = paf.t_len;
    flipped.q_st = paf.t_st;
    flipped.q_en = paf.t_en;
    // flip cigar
    let mut new_cigar = Vec::new();
    for opt in paf.cigar.into_iter() {
        let new_opt = match opt {
            Ins(l) => Del(*l),
            Del(l) => Ins(*l),
            _ => *opt,
        };
        new_cigar.push(new_opt);
    }
    if paf.strand == '-' {
        new_cigar.reverse();
    }
    flipped.cigar = CigarString(new_cigar);
    flipped
}

pub fn trim_paf_rec_to_rgn(rgn: &bed::Region, paf: &PafRecord) -> PafRecord {
    //let bytes = paf.cigar.as_bytes();
    //let cigar = CigarString::try_from(bytes).expect("Unable to parse cigar");

    // initalize a trimmed paf record
    let mut trimmed_paf = (*paf).clone();
    // check if we can return right away
    if paf.t_st >= rgn.st && paf.t_en <= rgn.en {
        //eprintln!("quick finish");
        return trimmed_paf;
    }
    // count until we get to the correct regions
    let mut t_pos = paf.t_st;
    let mut q_offset = 0;
    let mut q_offset_st = 0;
    let mut q_offset_en = 0;
    let mut new_cigar = Vec::new();
    for opt in paf.cigar.into_iter() {
        let moves_t = consumes_reference(&opt);
        let moves_q = consumes_query(&opt);
        let opt_len = opt.len() as u64;
        let mut new_opt_len = 0;
        let mut done = false;
        // incrment the ref and or query
        for _i in 0..opt_len {
            if moves_t {
                t_pos += 1;
            }
            if moves_q {
                q_offset += 1;
            }
            // add to the new cigar opt length
            if t_pos > rgn.st && t_pos <= rgn.en {
                new_opt_len += 1;
            }
            // record the starting positions
            if moves_t && t_pos == rgn.st {
                trimmed_paf.t_st = rgn.st;
                q_offset_st = q_offset;
            }
            // record the ending positions
            if moves_t && (t_pos == rgn.en || t_pos == paf.t_en) {
                trimmed_paf.t_en = t_pos;
                q_offset_en = q_offset;
                done = true;
                break;
            }
        }
        // add to the new cigar string.
        if new_opt_len != 0 {
            new_cigar.push(match opt {
                Match(_) => Match(new_opt_len),
                Ins(_) => Ins(new_opt_len),
                Del(_) => Del(new_opt_len),
                RefSkip(_) => RefSkip(new_opt_len),
                HardClip(_) => HardClip(new_opt_len),
                SoftClip(_) => SoftClip(new_opt_len),
                Pad(_) => Pad(new_opt_len),
                Equal(_) => Equal(new_opt_len),
                Diff(_) => Diff(new_opt_len),
            })
        }
        if done {
            break;
        }

        /*
        // record the starting positions
        if moves_t && t_pos >= rgn.st && t_pos - opt_len < rgn.st {
            trimmed_paf.t_st = rgn.st;
            q_offset_st = q_offset;
            if moves_q {
                q_offset_st -= t_pos - rgn.st;
            }
            eprintln!("set start\t{}\t{}", rgn.st, q_offset_st);
        } else if t_pos > rgn.st { // add to the cigar string
            new_cigar.push(opt.clone());
        }
        // record the ending positions
        if (moves_t && t_pos >= rgn.en && t_pos - opt_len < rgn.en) ||
            t_pos == paf.t_en {
            trimmed_paf.t_en = rgn.en;
            q_offset_en = q_offset;
            if moves_q {
                q_offset_en -= t_pos - rgn.en;
            }
            eprintln!("set end\t{}\t{}", rgn.en, q_offset_en);
            break;
        }*/
    }
    // fix strand if needed
    if paf.strand == '+' {
        trimmed_paf.q_st = paf.q_st + q_offset_st;
        trimmed_paf.q_en = paf.q_st + q_offset_en;
    } else {
        trimmed_paf.q_st = paf.q_en - q_offset_en + 1;
        trimmed_paf.q_en = paf.q_en - q_offset_st;
    }
    trimmed_paf.cigar = CigarString(new_cigar); // update to trimmed cigar.
    trimmed_paf
}

pub fn trim_paf_to_rgn(rgn: &bed::Region, paf: &[PafRecord], invert_query: bool) -> Vec<PafRecord> {
    let mut trimmed_paf = Vec::new();
    for rec in paf {
        if invert_query {
            //eprintln!("Inverting the PAF.");
            let rec = &paf_swap_query_and_target(rec);
            if paf_overlaps_target(rec, rgn) {
                trimmed_paf.push(trim_paf_rec_to_rgn(rgn, rec));
            }
        } else if paf_overlaps_target(rec, rgn) {
            trimmed_paf.push(trim_paf_rec_to_rgn(rgn, rec));
        }
    }
    trimmed_paf
}

/// # Example
/// ```
/// use rustybam::paf;
///
/// paf::read_paf_line("A 1 2 3 + B 1 2 3 10 11 60").unwrap();
/// let rec = paf::make_fake_paf_rec();
/// assert_eq!("4M1I1D3=", rec.cigar.to_string());
///
/// ```
pub fn read_paf_line(line: &str) -> PafResult<PafRecord> {
    let t: Vec<&str> = line.split_ascii_whitespace().collect();
    assert!(t.len() >= 12); // must have all required columns
                            // collect all the tags if any
    let mut tags = "".to_string();
    // find the cigar if it is there
    let mut cigar = "".to_string();
    let pattern = Regex::new("(..):(.):(.*)").unwrap();
    for token in t.iter().skip(12) {
        assert!(pattern.is_match(token));
        let caps = pattern.captures(token).unwrap();
        let tag = &caps[1];
        let value = &caps[3];
        if tag == "cg" {
            cigar = value.to_string();
        } else {
            tags.push('\t');
            tags.push_str(token);
        }
    }
    let cigar = cigar_from_str(cigar.as_str())?;
    // make the record
    let rec = PafRecord {
        q_name: t[0].to_string(),
        q_len: t[1].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        q_st: t[2].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        q_en: t[3].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        strand: t[4].parse::<char>().map_err(|_| Error::ParsePafColumn {})?,
        t_name: t[5].to_string(),
        t_len: t[6].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        t_st: t[7].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        t_en: t[8].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        nmatch: t[9].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        aln_len: t[10].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        mapq: t[11].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
        cigar,
        tags,
    };
    Ok(rec)
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
    for (index, line) in paf_file.lines().enumerate() {
        eprint!("\rReading PAF line: {}", index + 1);
        match read_paf_line(&line.unwrap()) {
            Ok(rec) => vec.push(rec),
            Err(_) => eprintln!("\nUnable to parse PAF record. Skipping line {}", index + 1),
        }
    }
    eprintln!();
    vec
}

pub fn make_fake_paf_rec() -> PafRecord {
    read_paf_line("Q 10 2 10 - T 20 12 20 3 9 60 cg:Z:4M1I1D3=").unwrap()
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
    #[test]
    /// this function tests that if we subset the ref and then invert that paf
    /// we can get back the original ref subset.
    fn check_inversable() {
        let paf = make_fake_paf_rec();
        let rgn = Region {
            name: "T".to_string(),
            st: 14,
            en: 17,
            id: "None".to_string(),
        };
        let trim_paf = trim_paf_rec_to_rgn(&rgn, &paf);
        eprintln!("{}", trim_paf);
        let q_rgn = Region {
            name: trim_paf.q_name.clone(),
            st: trim_paf.q_st,
            en: trim_paf.q_en,
            id: "None".to_string(),
        };
        let ref_paf = paf_swap_query_and_target(&trim_paf_rec_to_rgn(
            &q_rgn,
            &paf_swap_query_and_target(&trim_paf),
        ));
        // make sure we recreated what we wanted
        assert_eq!(ref_paf.t_st, trim_paf.t_st);
        assert_eq!(ref_paf.t_en, trim_paf.t_en);
        assert_eq!(ref_paf.q_st, trim_paf.q_st);
        assert_eq!(ref_paf.q_en, trim_paf.q_en);
        assert_eq!(ref_paf.cigar, trim_paf.cigar);
        eprintln!("{}", ref_paf);
    }
}