use super::bed;
//use rayon::prelude::*;
use num_format::{Locale, ToFormattedString};
use std::collections::HashMap;

pub fn bed_stats(bed: &str, readable: bool, column: Option<u8>) {
    let rgns = bed::parse_bed(bed);
    match column {
        Some(c) => {
            let mut dict = HashMap::new();
            for rgn in rgns {
                let (bp, n) = dict
                    .entry(rgn.get_column(c).clone())
                    .or_insert((0_u64, 0_u64));
                *bp += rgn.en - rgn.st;
                *n += 1;
            }
            log::trace!("{:?}", dict);
            for (key, (count, n)) in dict.iter_mut() {
                if readable {
                    println!(
                        "{}\t{}\t{}",
                        key,
                        count.to_formatted_string(&Locale::en),
                        n.to_formatted_string(&Locale::en)
                    )
                } else {
                    println!("{}\t{}\t{}", key, count, n);
                }
            }
        }
        None => {
            let n = rgns.len();
            let count: u64 = rgns.into_iter().map(|rgn| rgn.en - rgn.st).sum();
            if readable {
                println!(
                    "{}\t{}",
                    count.to_formatted_string(&Locale::en),
                    n.to_formatted_string(&Locale::en)
                )
            } else {
                println!("{}\t{}", count, n);
            }
        }
    }
}
