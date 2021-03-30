use bio::data_structures::suffix_array::{lcp, shortest_unique_substrings, suffix_array};

/// make an array of the shortest shared substring at each position in the text
/// # Example
/// ```
/// use rustybam::suns::*;
/// let text = b"GCTGCTA$";
/// let sus = get_shortest_subseq_size(text);
/// assert_eq!(
///        sus,
///        [
///            Some(4),
///            Some(3),
///            Some(2),
///            Some(4),
///            Some(3),
///            Some(2),
///            Some(1),
///            Some(1)
///        ]
///    );
/// eprintln!("testing");
/// ```
pub fn get_shortest_subseq_size(text: &[u8]) -> Vec<Option<usize>> {
    let pos = suffix_array(text);
    // obtain compressed LCP array
    let lcp = lcp(text, &pos);

    // calculate shortest unique substrings
    shortest_unique_substrings(&pos, &lcp)
}

/// Returns a vector of start and end coordiantes that are 100% made of SUNs
/// Coordinates are bed style [).
/// # Example
/// ```
/// use rustybam::suns::*;
/// let text = b"GCTGCTA$";
/// let sus = get_shortest_subseq_size(text);
/// for (start, end) in find_intervals(sus, 2){
///     println!("{}, {}", start, end);  
/// }
/// ```
pub fn find_intervals(sus: Vec<Option<usize>>, kmer_size: usize) -> Vec<(usize, usize)> {
    let mut vec = Vec::new();
    let mut i: usize = 0;
    while i < sus.len() {
        let start = i; // the start of this interval
        let mut cur: usize = sus[i].unwrap();
        // iterate until we do not find a unique substring
        while cur <= kmer_size && (i + 1) < sus.len() {
            i += 1;
            cur = sus[i].unwrap();
        }
        // add the interval to the return (+1 for bed interval)
        let end = i + 1;
        if end - 1 > start {
            vec.push((start, end));
        }
        // increment to next start
        i += 1;
    }
    vec
}

/// prints the interval that is a SUN in a pretty way.
pub fn pretty_interval_print(intervals: &[(usize, usize)], text: &[u8]) {
    println!("\n{}", "-".repeat(50));
    for (start, end) in intervals {
        println!("start:{}, end:{}", start, end);
        println!("{}", std::str::from_utf8(text).unwrap());
        for i in 0..text.len() {
            if i >= *start && i < *end {
                print!(".");
            } else {
                print!(" ");
            }
        }
        println!();
    }
    println!("{}\n", "-".repeat(50));
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_get_shortest_subseq_size() {
        let text = b"GCTGCTA$";
        let sus = get_shortest_subseq_size(text);
        let intervals = find_intervals(sus, 2);
        assert_eq!(intervals, [(2, 4), (5, 8)]);
        pretty_interval_print(&intervals, text);
    }
    #[test]
    fn test_get_shortest_subseq_size_two() {
        let text = b"AAAGGG$";
        let sus = get_shortest_subseq_size(text);
        let intervals = find_intervals(sus, 2);
        pretty_interval_print(&intervals, text);
    }
}
