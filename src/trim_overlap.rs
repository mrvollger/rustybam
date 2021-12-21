use super::paf::*;
use log;
use rust_htslib::bam::record::Cigar::*;
use std::cmp::{max, min};

fn score_of_qpos(rec: &PafRecord, pos: u64) -> i32 {
    let aln_idx = rec.qpos_to_idx(pos).unwrap();
    match rec.long_cigar[aln_idx] {
        Equal(_) => 1,
        Ins(_) | Del(_) => -1,
        _ => -1,
    }
}

/// Example
/// ```
/// use rustybam::trim_overlap::*;
/// use rustybam::paf::*;
///
/// let mut left = PafRecord::new("Q 10 0 10 + T 20 0 10 3 9 60 cg:Z:7=1X2=").unwrap();
/// left.aligned_pairs();
/// let mut right = PafRecord::new("Q 10 5 10 - T 20 10 15 3 9 60 cg:Z:3=1X1=").unwrap();
/// right.aligned_pairs();
///
/// trim_overlapping_pafs(&mut left, &mut right);
///
/// assert_eq!(left.cigar.to_string(), "7=");
/// assert_eq!(right.cigar.to_string(), "3=");
/// ```
pub fn trim_overlapping_pafs(left: &mut PafRecord, right: &mut PafRecord) {
    let st_ovl = max(left.q_st, right.q_st);
    let en_ovl = min(left.q_en, right.q_en);

    let mut l_score = vec![0];
    let mut r_score = vec![];
    for pos in st_ovl..en_ovl {
        eprintln!("{}", pos);
        let l_s = score_of_qpos(left, pos);
        let r_s = score_of_qpos(right, pos);
        l_score.push(l_s);
        r_score.push(r_s);
    }
    r_score.push(0); //
                     // make cumulative from left
    l_score.iter_mut().fold(0, |acc, x| {
        *x += acc;
        *x
    });
    // make cumulative from right
    r_score.iter_mut().rev().fold(0, |acc, x| {
        *x += acc;
        *x
    });
    eprintln!("{:?}\n{:?}", l_score, r_score);

    // determine the split point
    let mut max_idx: u64 = 0;
    let mut max = 0;
    for (idx, (l, r)) in l_score.iter().zip(r_score.iter()).enumerate() {
        if l + r > max {
            max = l + r;
            max_idx = idx as u64;
        }
    }
    left.truncate_record_by_query(left.q_st, st_ovl + max_idx);
    right.truncate_record_by_query(st_ovl + max_idx, right.q_en);

    log::info!(
        "Split position was found to be {} bases after the overlap start ({}).",
        max_idx,
        st_ovl
    );
}
