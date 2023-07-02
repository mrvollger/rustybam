use super::paf;
use bio_types::strand::ReqStrand::*;
use colored::Colorize;
use rust_htslib::bam::record::{Cigar::*, CigarStringView};
use rust_htslib::bam::Header;
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::Record;
use std::convert::TryFrom;
use std::fmt;
use std::str;
#[derive(Default)]
pub struct Stats {
    pub q_nm: String,
    pub q_len: i64,
    pub q_st: i64,
    pub q_en: i64,
    pub r_nm: String,
    pub r_len: i64,
    pub r_st: i64,
    pub r_en: i64,
    pub strand: char,
    pub equal: u32,
    pub diff: u32,
    pub ins: u32,
    pub del: u32,
    pub matches: u32,
    pub ins_events: u32,
    pub del_events: u32,
    pub id_by_all: f32,
    pub id_by_events: f32,
    pub id_by_matches: f32,
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} {} {} {} {} {} {}",
            self.r_nm, self.r_st, self.r_en, self.strand, self.q_nm, self.q_st, self.q_en
        )
    }
}

pub fn stats_from_paf(paf: paf::PafRecord) -> Stats {
    //let paf = paf::read_paf_line(line).unwrap();
    let mut stats = Stats::default();
    add_stats_from_cigar(&CigarStringView::new(paf.cigar, 0), &mut stats);
    stats.r_nm = paf.t_name;
    stats.r_len = paf.t_len as i64;
    stats.r_st = paf.t_st as i64;
    stats.r_en = paf.t_en as i64;
    stats.q_nm = paf.q_name;
    stats.q_len = paf.q_len as i64;
    stats.q_st = paf.q_st as i64;
    stats.q_en = paf.q_en as i64;
    stats.strand = paf.strand;
    stats
}

pub fn add_stats_from_cigar(cigar: &CigarStringView, stats: &mut Stats) {
    // iterate over cigar
    for opt in cigar {
        match opt {
            Del(val) => {
                stats.del_events += 1;
                stats.del += val
            }
            Ins(val) => {
                stats.ins_events += 1;
                stats.ins += val
            }
            Equal(val) => stats.equal += val,
            Diff(val) => stats.diff += val,
            Match(val) => {
                stats.diff += val;
                stats.matches += val
            }
            _ => (),
        }
    }

    // make the summary stats
    stats.id_by_all =
        100.0 * stats.equal as f32 / (stats.equal + stats.diff + stats.del + stats.ins) as f32;
    stats.id_by_events = 100.0 * stats.equal as f32
        / (stats.equal + stats.diff + stats.del_events + stats.ins_events) as f32;
    stats.id_by_matches = 100.0 * stats.equal as f32 / (stats.equal + stats.diff) as f32;

    // print warnings if the cigar does not use =/X
    if stats.matches > 0 {
        eprint!(
            "\r{} {} {}{}",
            "\u{26A0} warning:".bold().yellow(),
            "cigar string contains".yellow(),
            "'M'".bold().red(),
            ", assuming mismatch.".yellow()
        );
    }
}

pub fn cigar_stats(rec: Record, header: &Header) -> Stats {
    let cigar = rec.cigar();
    let bam_head = HeaderView::from_header(header);
    let r_nm = str::from_utf8(bam_head.tid2name(rec.tid() as u32)).unwrap();
    let r_len = bam_head.target_len(rec.tid() as u32).unwrap();

    // initalize output information
    let mut stats = Stats {
        r_nm: r_nm.to_string(),
        r_len: r_len as i64,
        r_st: rec.pos(),
        r_en: cigar.end_pos(),
        q_nm: std::str::from_utf8(rec.qname()).unwrap().to_string(),
        q_len: 0,
        q_st: 0,
        q_en: 0,
        strand: rec
            .strand()
            .strand_symbol()
            .chars()
            .next()
            .expect("string is empty"),
        equal: 0,
        diff: 0,
        ins: 0,
        del: 0,
        ins_events: 0,
        del_events: 0,
        matches: 0,
        id_by_matches: 0.0,
        id_by_events: 0.0,
        id_by_all: 0.0,
    };

    // get the query coordinates
    stats.q_st = cigar.leading_hardclips() + cigar.leading_softclips();
    stats.q_en = cigar.leading_hardclips()
        + 1
        + cigar
            .read_pos(stats.r_en as u32 - 1, false, false)
            .unwrap()
            .unwrap() as i64; // the plus 1 is to make the end exclusive [x, y) format

    stats.q_len = cigar.leading_hardclips()
        + i64::try_from(rec.seq_len()).unwrap()
        + cigar.trailing_hardclips();
    // fix query coordinates if rc
    if rec.strand() == Reverse {
        let temp = stats.q_st;
        stats.q_st = stats.q_len - stats.q_en;
        stats.q_en = stats.q_len - temp;
    }

    // read the cigar string and add in the stats
    add_stats_from_cigar(&cigar, &mut stats);

    //eprintln!("{} {} {} {}", stats.equal, stats.diff, stats.ins, stats.del);
    // println!("{}", stats);
    stats
}

/// print stats
pub fn print_cigar_stats_header(qbed: bool) {
    if qbed {
        print!("#query_name\tquery_start\tquery_end\tquery_length\t");
        print!("strand\t");
        print!("reference_name\treference_start\treference_end\treference_length\t");
    } else {
        print!("#reference_name\treference_start\treference_end\treference_length\t");
        print!("strand\t");
        print!("query_name\tquery_start\tquery_end\tquery_length\t");
    }
    println!("perID_by_matches\tperID_by_events\tperID_by_all\tmatches\tmismatches\tdeletion_events\tinsertion_events\tdeletions\tinsertions");
}

/// print cigar stats from a bam
pub fn print_cigar_stats(stats: Stats, qbed: bool) {
    if qbed {
        print!(
            "{}\t{}\t{}\t{}\t",
            stats.q_nm, stats.q_st, stats.q_en, stats.q_len
        );
        print!("{}\t", stats.strand);
        print!(
            "{}\t{}\t{}\t{}\t",
            stats.r_nm, stats.r_st, stats.r_en, stats.r_len
        );
    } else {
        print!(
            "{}\t{}\t{}\t{}\t",
            stats.r_nm, stats.r_st, stats.r_en, stats.r_len
        );
        print!("{}\t", stats.strand);
        print!(
            "{}\t{}\t{}\t{}\t",
            stats.q_nm, stats.q_st, stats.q_en, stats.q_len
        );
    }

    print!(
        "{}\t{}\t{}\t",
        stats.id_by_matches, stats.id_by_events, stats.id_by_all
    );
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        stats.equal, stats.diff, stats.del_events, stats.ins_events, stats.del, stats.ins
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_stats_from_test_file() {
        use rust_htslib::{bam, bam::Read};
        let mut bam = bam::Reader::from_path(".test/asm_small.bam").unwrap();
        let bam_header = bam::Header::from_template(bam.header());
        bam.set_threads(4).unwrap();
        for fetch in bam.records() {
            let stats = cigar_stats(fetch.unwrap(), &bam_header);
            print_cigar_stats(stats, false);
        }
    }
    #[test]
    fn test_add_cigar_stats() {
        use rust_htslib::bam::record::CigarString;
        let cigar = CigarString::try_from("10=10X").unwrap();
        let view = CigarStringView::new(cigar, 0);
        let mut stats = Stats::default();
        add_stats_from_cigar(&view, &mut stats);
        assert!((50.0 - stats.id_by_all).abs() < 1e-10);
    }
}
