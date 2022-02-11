use super::paf::*;
use log;
use rust_htslib::bam::record::Cigar::*;
use std::cmp::{max, min};

fn score_of_qpos(
    rec: &PafRecord,
    pos: u64,
    match_score: i32,
    diff_score: i32,
    indel_score: i32,
) -> i32 {
    let aln_idx = rec.qpos_to_idx(pos).unwrap();
    match rec.long_cigar[aln_idx] {
        Equal(_) => match_score,
        Ins(_) | Del(_) => -indel_score,
        _ => -diff_score,
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
/// trim_overlapping_pafs(&mut left, &mut right, 1 ,1 ,1);
///
/// assert_eq!(left.cigar.to_string(), "7=");
/// assert_eq!(right.cigar.to_string(), "3=");
/// ```
pub fn trim_overlapping_pafs(
    left: &mut PafRecord,
    right: &mut PafRecord,
    match_score: i32,
    diff_score: i32,
    indel_score: i32,
) {
    let st_ovl = max(left.q_st, right.q_st);
    let en_ovl = min(left.q_en, right.q_en);
    log::info!("Number of overlapping bases {}", en_ovl - st_ovl);

    let mut l_score = vec![0];
    let mut r_score = vec![];

    for pos in st_ovl..en_ovl {
        let l_s = score_of_qpos(left, pos, match_score, diff_score, indel_score);
        let r_s = score_of_qpos(right, pos, match_score, diff_score, indel_score);
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
        "Split position was found to be {} bases after the overlap start ({}) with a score of {}.",
        max_idx,
        st_ovl,
        max
    );
}

/// Example
/// ```
/// use rustybam::trim_overlap::*;
/// use rustybam::paf::*;
///
/// ```

#[cfg(test)]
mod tests {
    use super::*;

    fn pretty_print_paf(paf: &Paf) {
        println!("\nQUERY");
        for rec in &paf.records {
            for _i in 0..rec.q_st {
                print!(" ");
            }
            if rec.strand == '-' {
                rec.long_cigar.iter().rev().for_each(|x| {
                    let opt = format!("{}", x).strip_prefix("1").unwrap().to_string();
                    print!("{}", opt);
                });
            } else {
                rec.long_cigar.iter().for_each(|x| {
                    let opt = format!("{}", x).strip_prefix("1").unwrap().to_string();
                    print!("{}", opt);
                });
            }

            for _i in 0..(rec.q_len - rec.q_en + 5) {
                print!(" ");
            }
            println!("{}-{}", rec.q_st, rec.q_en);
        }

        println!("REFERENCE");
        for rec in &paf.records {
            for _i in 0..rec.t_st {
                print!(" ");
            }
            rec.long_cigar.iter().for_each(|x| {
                let opt = format!("{}", x).strip_prefix("1").unwrap().to_string();
                print!("{}", opt);
            });

            for _i in 0..(rec.t_len - rec.t_en + 5) {
                print!(" ");
            }
            println!("{}-{}", rec.t_st, rec.t_en);
        }
    }

    #[test]
    fn test_inversion_trimming() {
        let mut left = PafRecord::new("Q 20 0 10 + T 20 0 10 3 9 60 cg:Z:7=1X2=").unwrap();
        left.aligned_pairs();
        left.check_integrity().unwrap();

        let mut center = PafRecord::new("Q 20 4 15 - T 20 5 16 3 9 60 cg:Z:3=1X3=1M1X2=").unwrap();
        center.aligned_pairs();
        center.check_integrity().unwrap();

        let mut right =
            PafRecord::new("Q 20 10 20 + T 20 10 20 3 9 60 cz:Z:10= cg:Z:2=2X2=2X2=").unwrap();
        right.aligned_pairs();
        right.check_integrity().unwrap();

        println!("{:?}", center.tpos_aln);
        println!("{:?}", center.qpos_aln);

        let mut paf = Paf::new();
        paf.records = vec![left, center, right];

        pretty_print_paf(&paf);
        paf.overlapping_paf_recs(1, 1, 1);
        pretty_print_paf(&paf);

        let expected_cigars = vec!["7=", "2=1X3=1M", "2=2X2="];
        for (i, rec) in paf.records.iter().enumerate() {
            assert_eq!(rec.cigar.to_string(), expected_cigars[i]);
        }
    }
}
