use super::bed;
//use rayon::prelude::*;
use std::collections::HashMap;

pub fn bed_stats(bed: &str, readable: bool, column: Option<u8>) {
    let rgns = bed::parse_bed(bed);
    match column {
        Some(c) => {
            let mut dict = HashMap::new();
            for rgn in rgns {
                let (bp, n) = dict
                    .entry(rgn.get_column(c).clone())
                    .or_insert((0_f32, 0_u64));
                *bp += (rgn.en - rgn.st) as f32;
                *n += 1;
            }
            log::trace!("{:?}", dict);
            for (key, (count, n)) in dict.iter_mut() {
                if readable {
                    *count /= 1e6;
                }
                println!("{}\t{}\t{}", key, count, n);
            }
        }
        None => {
            let n = rgns.len();
            let mut count: f32 = rgns.into_iter().map(|rgn| (rgn.en - rgn.st) as f32).sum();
            if readable {
                count /= 1e6;
            }
            println!("{}\t{}", count, n);
        }
    }
}
