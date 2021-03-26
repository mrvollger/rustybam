use std::str;
use rust_htslib::bam::record::Cigar::*;
use rust_htslib::bam::Record;
use bio_types::strand::ReqStrand;
use bio_types::strand::ReqStrand::*;
use std::fmt;
use colored::Colorize;
use std::convert::TryFrom;
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::Header;


pub struct Stats{
  pub tid : i32,
  pub q_nm : String,
  pub q_len : i64,
  pub q_st : i64,
  pub q_en : i64,
  pub r_st : i64,
  pub r_en : i64,
  pub strand : ReqStrand,
  pub equal : u32,
  pub diff : u32,
  pub ins : u32,
  pub del : u32,
  pub matches : u32,
  pub ins_events : u32,
  pub del_events : u32,
  pub id_by_all : f32,
  pub id_by_events : f32,
  pub id_by_matches : f32,
}

impl fmt::Display for Stats {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f,
          "{} {} {} {} {} {} {}", 
          self.tid, self.r_st, self.r_en,
          self.strand.strand_symbol(), 
          self.q_nm, self.q_st, self.q_en
        )
  }
}

pub fn cigar_stats(mut rec : Record) -> Stats{
  let cigar = rec.cigar();

  // initalize output information 
  let mut stats = Stats{
    tid : rec.tid(), r_st : rec.pos(), r_en : cigar.end_pos(),
    q_nm : std::str::from_utf8(rec.qname()).unwrap().to_string(), q_len : 0,    
    q_st : 0, q_en : 0, strand : rec.strand(),
    equal : 0, diff : 0,  ins : 0, del : 0, ins_events : 0, del_events : 0, matches : 0,
    id_by_matches : 0.0, id_by_events : 0.0, id_by_all : 0.0 
  };

  // get the query coordinates
  stats.q_st = cigar.leading_hardclips() + cigar.leading_softclips(); 
  stats.q_en = cigar.leading_hardclips() + 1 + // the plus 1 is to make the end exclusive [x, y) format 
                cigar.read_pos(stats.r_en as u32 - 1, false, false).unwrap().unwrap() as i64;
  stats.q_len = cigar.leading_hardclips() + i64::try_from(rec.seq_len()).unwrap() + cigar.trailing_hardclips();
  
  // fix query coordinates if rc
  match rec.strand(){
    Reverse => {
      let temp = stats.q_st;
      stats.q_st = stats.q_len - stats.q_en;
      stats.q_en = stats.q_len - temp;
    },
    _ => (),
  }

  // count up the cigar operations 
  for opt in &cigar  {
    match opt {
      Del(val) => {stats.del_events += 1; stats.del += val},
      Ins(val) => {stats.ins_events += 1; stats.ins += val},
      Equal(val) => {stats.equal += val},
      Diff(val) => {stats.diff += val},
      Match(val) => {stats.diff += val; stats.matches += val},
      _ => (),
    }
  }

  stats.id_by_all = 100.0 * stats.equal as f32 / (stats.equal + stats.diff + stats.del + stats.ins) as f32;
  stats.id_by_events = 100.0 * stats.equal as f32 / (stats.equal + stats.diff + stats.del_events + stats.ins_events) as f32;
  stats.id_by_matches = 100.0 * stats.equal as f32 / (stats.equal + stats.diff) as f32;

  if stats.matches > 0 {
    eprint!("\r{} {} {}{}",
      "\u{26A0} warning:".bold().yellow(), 
      "cigar string contains".yellow(),
      "'M'".bold().red(),
      ", assuming mismatch.".yellow());
  }
  //eprintln!("{} {} {} {}", stats.equal, stats.diff, stats.ins, stats.del);
  // println!("{}", stats);
  stats
}  

/// print stats
pub fn print_cigar_stats_header(qbed : bool){
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
pub fn print_cigar_stats( stats : Stats, qbed : bool, header : &Header ){
  let bam_head = HeaderView::from_header(header);
  let r_nm = str::from_utf8( bam_head.tid2name(stats.tid as u32) ).unwrap();
  let r_len = bam_head.target_len(stats.tid as u32).unwrap();

  if qbed {
    print!("{}\t{}\t{}\t{}\t", stats.q_nm, stats.q_st, stats.q_en, stats.q_len);
    print!("{}\t", stats.strand.strand_symbol());
    print!("{}\t{}\t{}\t{}\t", r_nm, stats.r_st, stats.r_en, r_len);
  } else {
    print!("{}\t{}\t{}\t{}\t", r_nm, stats.r_st, stats.r_en, r_len);
    print!("{}\t", stats.strand.strand_symbol());
    print!("{}\t{}\t{}\t{}\t", stats.q_nm, stats.q_st, stats.q_en, stats.q_len);
  }
  
  print!("{}\t{}\t{}\t", stats.id_by_matches, stats.id_by_events, stats.id_by_all);
  print!("{}\t{}\t{}\t{}\t{}\t{}\n", stats.equal, stats.diff, stats.del_events, stats.ins_events,  stats.del, stats.ins);
}

#[cfg(test)]
mod tests {
  use super::*;

    #[test]
    fn test_cigar_stats_from_test_file() {
      use rust_htslib::{bam, bam::Read};
      let mut bam = bam::Reader::from_path("test/asm_small.bam").unwrap();
      let bam_header = bam::Header::from_template(bam.header());
      bam.set_threads(4).unwrap();
      for fetch in bam.records(){
        let stats = cigar_stats(fetch.unwrap()); 
        print_cigar_stats(stats, false, &bam_header);
      }
    }
}
