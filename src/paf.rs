use super::bed;
use core::{fmt, panic};
use regex::Regex;
use rust_htslib::bam::record::Cigar::*;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::record::*;
use std::collections::HashMap;
use std::fs;
use std::io;
use std::io::BufRead;
use std::str::FromStr;
use std::usize;

#[derive(Debug)]
pub enum Error {
    PafParseCigar { msg: String },
    PafParseCS { msg: String },
    ParseIntError { msg: String },
    ParsePafColumn {},
}
type PafResult<T> = Result<T, crate::paf::Error>;

#[derive(Debug)]
pub struct Paf<'a> {
    pub records: Vec<PafRecord>,
    pub records_by_contig: HashMap<String, Vec<&'a PafRecord>>,
}

impl<'a> Paf<'a> {
    fn new() -> Paf<'a> {
        Paf {
            records: Vec::new(),
            records_by_contig: HashMap::new(),
        }
    }
    /// read in the paf from a file pass "-" for stdin
    /// # Example
    /// ```
    /// use rustybam::paf;
    /// use std::fs::File;
    /// use std::io::*;
    /// let mut paf = paf::Paf::from_file(".test/asm_small.paf");
    /// assert_eq!(paf.records.len(), 249);
    ///
    /// ```
    pub fn from_file(file_name: &str) -> Paf<'a> {
        // open the paf file
        let paf_file: Box<dyn io::Read> = match file_name {
            "-" => Box::new(io::stdin()),
            _ => Box::new(fs::File::open(file_name).expect("Unable to open paf file")),
        };
        let mut paf = Paf::new();
        let mut records = Vec::new();
        // read the paf recs into a vector
        for (index, line) in io::BufReader::new(paf_file).lines().enumerate() {
            match PafRecord::new(&line.unwrap()) {
                Ok(rec) => {
                    eprint!("\rReading PAF entry # {}", index);
                    records.push(rec)
                }
                Err(_) => eprintln!("\nUnable to parse PAF record. Skipping line {}", index + 1),
            }
        }
        eprintln!();
        paf.records = records;
        paf
    }
    /*
    pub fn populate_hash(&'a mut self) {
        // read pafs into a hashtable based on name
        for rec in &self.records[..] {
            if self.records_by_contig.contains_key(&rec.t_name) {
                self.records_by_contig
                    .insert(rec.t_name.clone(), Vec::new());
            }
            self.records_by_contig
                .get_mut(&rec.t_name)
                .unwrap()
                .push(&rec);
        }
    }

    // calculate infor for different records
    pub fn make_position_index(&mut self) {
        eprintln!("Making paf index:");
        self.records
            .par_iter_mut()
            .for_each(|rec| rec.aligned_pairs());
    }
    */
}

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
    pub tpos_aln: Vec<u64>,
    pub qpos_aln: Vec<u64>,
    pub long_cigar: CigarString,
    pub id: String,
}

impl PafRecord {
    /// # Example
    /// ```
    /// use rustybam::paf;
    /// let _paf = paf::PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60").unwrap();
    /// let rec = paf::make_fake_paf_rec();
    /// assert_eq!("4M1I1D3=", rec.cigar.to_string());
    ///
    /// ```
    pub fn new(line: &str) -> PafResult<PafRecord> {
        let t: Vec<&str> = line.split_ascii_whitespace().collect();
        assert!(t.len() >= 12); // must have all required columns
                                // collect all the tags if any
        let mut tags = "".to_string();
        // find the cigar if it is there
        let mut cigar = CigarString(vec![]);
        let pattern = Regex::new("(..):(.):(.*)").unwrap();
        for token in t.iter().skip(12) {
            assert!(pattern.is_match(token));
            let caps = pattern.captures(token).unwrap();
            let tag = &caps[1];
            let value = &caps[3];
            // TODO fix cs string parsing when both cigar and cs are there.
            if tag == "cs" {
                cigar = cs_to_cigar(value)?;
                eprintln!("cs parsed");
            } else if tag == "cg" && cigar.len() == 0 {
                cigar = cigar_from_str(value)?;
            } else {
                tags.push('\t');
                tags.push_str(token);
            }
        }

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
            tpos_aln: Vec::new(),
            qpos_aln: Vec::new(),
            long_cigar: CigarString(Vec::new()),
            id: "".to_string(),
        };
        Ok(rec)
    }

    pub fn small_copy(&self) -> PafRecord {
        PafRecord {
            q_name: self.q_name.clone(),
            q_len: self.q_len,
            q_st: self.q_st,
            q_en: self.q_en,
            strand: self.strand,
            t_name: self.t_name.clone(),
            t_len: self.t_len,
            t_st: self.t_st,
            t_en: self.t_en,
            nmatch: self.nmatch,
            aln_len: self.aln_len,
            mapq: self.mapq,
            cigar: CigarString(Vec::new()),
            tags: self.tags.clone(),
            tpos_aln: Vec::new(),
            qpos_aln: Vec::new(),
            long_cigar: CigarString(Vec::new()),
            id: self.id.clone(),
        }
    }

    /// This function adds matching alignment positions and cigar operations
    pub fn aligned_pairs(&mut self) {
        let mut t_pos = self.t_st as i64 - 1;
        let mut q_pos = self.q_st as i64 - 1;
        let mut long_cigar = Vec::new();
        if self.strand == '-' {
            q_pos = self.q_en as i64; // ends are not inclusive
        }

        for opt in self.cigar.into_iter() {
            let moves_t = consumes_reference(opt);
            let moves_q = consumes_query(opt);
            let opt_len = opt.len() as u64;
            // incrment the ref and or query
            for _i in 0..opt_len {
                long_cigar.push(update_cigar_opt_len(opt, 1));
                if moves_t {
                    t_pos += 1;
                }
                if moves_q && self.strand == '+' {
                    q_pos += 1;
                }
                if moves_q && self.strand == '-' {
                    q_pos -= 1;
                }
                self.tpos_aln.push(t_pos as u64);
                self.qpos_aln.push(q_pos as u64);
            }
        }

        self.long_cigar = CigarString(long_cigar);
    }

    // Does a binary search on the alignment to find its index in the alignment
    pub fn tpos_to_idx(&self, tpos: u64) -> usize {
        /*eprintln!(
            "{}<{}<{},  {}-{}",
            self.t_st,
            tpos,
            self.t_en,
            self.tpos_aln[0],
            self.tpos_aln[self.tpos_aln.len() - 1]
        );*/
        let mut idx = self
            .tpos_aln
            .binary_search(&tpos)
            .expect("must be within target alignment positions");
        while idx + 1 < self.tpos_aln.len() && self.tpos_aln[idx + 1] == tpos {
            idx += 1;
        }
        idx
    }

    // subset cigar on coordiantes
    pub fn subset_cigar(&self, start_idx: usize, end_idx: usize) -> CigarString {
        // update the cigar string
        let mut new_cigar = vec![];
        for opt in &self.long_cigar.0[start_idx..end_idx + 1] {
            new_cigar.push(*opt);
        }
        CigarString(new_cigar)
    }

    pub fn collapse_long_cigar(cigar: &CigarString) -> CigarString {
        let mut rtn = Vec::new();
        let mut pre_opt = cigar.0[0];
        let mut pre_len = 1;
        let mut idx = 1;
        while idx < cigar.len() {
            let cur_opt = cigar.0[idx];
            if std::mem::discriminant(&cur_opt) == std::mem::discriminant(&pre_opt) {
                pre_len += 1;
            } else {
                rtn.push(update_cigar_opt_len(&pre_opt, pre_len));
                pre_opt = cur_opt;
                pre_len = 1;
            }
            idx += 1;
        }
        rtn.push(update_cigar_opt_len(&pre_opt, pre_len));
        CigarString(rtn)
    }

    pub fn paf_overlaps_rgn(&self, rgn: &bed::Region) -> bool {
        if self.t_name != rgn.name {
            return false;
        }
        self.t_en > rgn.st && self.t_st < rgn.en
    }

    pub fn remove_trailing_indels(&mut self) {
        let st_opt = *self.cigar.first().unwrap();
        let en_opt = *self.cigar.last().unwrap();

        let mut remove_st = 0;
        if matches!(st_opt, Ins(_) | Del(_)) {
            remove_st = st_opt.len();
            self.cigar = CigarString(self.cigar.0[1..].to_vec());
        }
        let mut remove_en = 0;
        if matches!(en_opt, Ins(_) | Del(_)) {
            remove_en = en_opt.len();
            //self.cigar = CigarString(self.cigar.0[0..self.cigar.len() - 1].to_vec());
            self.cigar.0.truncate(self.cigar.len() - 1);
        }

        if matches!(st_opt, Del(_)) || matches!(en_opt, Del(_)) {
            self.t_st += remove_st as u64;
            self.t_en -= remove_en as u64;
        }
        // change the qpos removal
        if self.strand == '-' {
            std::mem::swap(&mut remove_st, &mut remove_en);
        }
        // fix the query positions that need to be
        if matches!(st_opt, Ins(_)) || matches!(en_opt, Ins(_)) {
            self.q_st += remove_st as u64;
            self.q_en -= remove_en as u64;
        }
    }
}

impl fmt::Display for PafRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tid:Z:{}\tcg:Z:{}",
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
            self.id,
            self.cigar.to_string(),
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
/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::Cigar::*;
/// assert_eq!(Diff(5), paf::update_cigar_opt_len(&Diff(10), 5));
/// assert_eq!(Diff(10), paf::update_cigar_opt_len(&Diff(1), 10));
/// ```
pub fn update_cigar_opt_len(opt: &Cigar, new_opt_len: u32) -> Cigar {
    match opt {
        Match(_) => Match(new_opt_len),
        Ins(_) => Ins(new_opt_len),
        Del(_) => Del(new_opt_len),
        RefSkip(_) => RefSkip(new_opt_len),
        HardClip(_) => HardClip(new_opt_len),
        SoftClip(_) => SoftClip(new_opt_len),
        Pad(_) => Pad(new_opt_len),
        Equal(_) => Equal(new_opt_len),
        Diff(_) => Diff(new_opt_len),
    }
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

/// Basically swaps the query and the reference in a cigar
pub fn cigar_swap_target_query(cigar: &CigarString, strand: char) -> CigarString {
    // flip cigar
    let mut new_cigar = Vec::new();
    for opt in cigar.into_iter() {
        let new_opt = match opt {
            Ins(l) => Del(*l),
            Del(l) => Ins(*l),
            _ => *opt,
        };
        new_cigar.push(new_opt);
    }
    if strand == '-' {
        new_cigar.reverse();
    }
    CigarString(new_cigar)
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

    // flip the index
    std::mem::swap(&mut flipped.qpos_aln, &mut flipped.tpos_aln);

    if paf.strand == '-' {
        flipped.qpos_aln.reverse();
        flipped.tpos_aln.reverse();
    }
    // flip the cigar
    flipped.cigar = cigar_swap_target_query(&paf.cigar, paf.strand);
    flipped.long_cigar = cigar_swap_target_query(&paf.long_cigar, paf.strand);
    flipped
}

pub fn old_trim_paf_rec_to_rgn(rgn: &bed::Region, paf: &PafRecord) -> PafRecord {
    // initalize a trimmed paf record
    let mut trimmed_paf = (*paf).clone();
    // check if we can return right away
    if paf.t_st >= rgn.st && paf.t_en <= rgn.en {
        return trimmed_paf;
    }
    // count until we get to the correct regions
    let mut t_pos = paf.t_st;
    let mut q_offset = 0;
    let mut q_offset_st = 0;
    let mut q_offset_en = 0;
    let mut new_cigar = Vec::new();
    trimmed_paf.nmatch = 0;
    trimmed_paf.aln_len = 0;
    for opt in paf.cigar.into_iter() {
        let moves_t = consumes_reference(opt);
        let moves_q = consumes_query(opt);
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
                trimmed_paf.aln_len += 1;
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
                Equal(_) => {
                    trimmed_paf.nmatch += new_opt_len as u64;
                    Equal(new_opt_len)
                }
                Diff(_) => Diff(new_opt_len),
            })
        }
        if done {
            break;
        }
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

pub fn make_fake_paf_rec() -> PafRecord {
    let mut rtn = PafRecord::new("Q 10 2 10 - T 20 12 20 3 9 60 cg:Z:4M1I1D3=").unwrap();
    rtn.aligned_pairs();
    rtn
}

/// # Example
/// ```
/// use rust_htslib::bam::record::Cigar::*;
/// use rustybam::paf;
/// let cigar = paf::cs_to_cigar(":10=ACGTN+acgtn-acgtn*at=A").unwrap();
/// assert_eq!(cigar[0], Equal(10));
/// assert_eq!(cigar[1], Equal(5));
/// assert_eq!(cigar[2], Ins(5));
/// assert_eq!(cigar[3], Del(5));
/// assert_eq!(cigar[4], Diff(1));
/// assert_eq!(cigar[5], Equal(1));
/// ```
pub fn cs_to_cigar(cs: &str) -> PafResult<CigarString> {
    let bytes = cs.as_bytes();
    let length = bytes.len();
    let mut i = 0;
    let mut cigar = vec![];
    while i < length {
        let cs_opt = bytes[i];
        let mut l = 0;
        // get past the opt and to the information
        i += 1;
        let opt = match cs_opt {
            b'=' => {
                while let b'A' | b'C' | b'G' | b'T' | b'N' = bytes[i] {
                    i += 1;
                    l += 1;
                    if i == length {
                        break;
                    }
                }
                Cigar::Equal(l)
            }
            b':' => {
                let mut j = i;
                while j < length && bytes[j].is_ascii_digit() {
                    j += 1;
                }
                l = u32::from_str(&cs[i..j]).map_err(|_| Error::ParseIntError {
                    msg: format!("Expected integer, got {}", cs[i..j].to_string()),
                })?;
                i += j - 1;
                Cigar::Equal(l)
            }
            b'*' => {
                i += 2;
                Cigar::Diff(1)
            }
            b'+' | b'-' => {
                while let b'a' | b'c' | b'g' | b't' | b'n' = bytes[i] {
                    i += 1;
                    l += 1;
                    if i == length {
                        break;
                    }
                }
                match cs_opt {
                    b'+' => Cigar::Ins(l),
                    b'-' => Cigar::Del(l),
                    _ => panic!("should be impossible + or - needed"),
                }
            }
            b'~' => {
                return Err(Error::PafParseCS {
                    msg: "Splice operations not yet supported.".to_string(),
                });
            }
            _ => {
                return Err(Error::PafParseCS {
                    msg: format!("Unexpected operator in the cs string: {}", cs_opt as char),
                });
            }
        };
        cigar.push(opt);
    }
    Ok(CigarString(cigar))
}
