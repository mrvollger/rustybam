use bio::alphabets::dna::revcomp;
use bio::data_structures::suffix_array::{lcp, shortest_unique_substrings, suffix_array};
use bio::io::fasta;
use std::{fs::File, usize};

static END_CHAR: u8 = b'$';
static END_CHAR_STR: &str = "$";

pub struct Genome {
    pub names: Vec<String>,
    pub starts: Vec<usize>,
    pub ends: Vec<usize>,
    pub seq: Vec<u8>,
    pub length: usize,
}

impl Genome {
    pub fn new(records: fasta::Records<File>) -> Genome {
        let mut genome = Genome {
            names: Vec::new(),
            starts: Vec::new(),
            ends: Vec::new(),
            seq: Vec::new(),
            length: 0,
        };
        let mut seq = Vec::new();
        for rec in records {
            let rec = rec.unwrap();
            // define the starts and names of contigs w.r.t the flat genome
            genome.starts.push(seq.len());
            genome.names.push(rec.id().to_string());
            // add the sequence
            seq.append(&mut rec.seq().to_ascii_uppercase());
            seq.push(END_CHAR);
            genome.ends.push(seq.len());
        }
        genome.length = seq.len();
        seq.append(&mut revcomp(&seq));
        seq.push(END_CHAR);
        genome.seq = seq;
        genome
    }

    /// # Example
    /// ```
    /// use rustybam::suns::*;
    /// let genome = Genome::from_file("test/test.fa");
    /// ```
    pub fn from_file(fastafile: &str) -> Genome {
        let records = fasta::Reader::from_file(fastafile)
            .expect("Unable to read")
            .records();
        Genome::new(records)
    }
    /// make an array of the shortest shared substring at each position in the text
    /// # Example
    /// ```
    /// use rustybam::suns::*;
    /// let text = b"GCTGCTA$";
    /// let sus = Genome::get_shortest_subseq_size(text);
    /// assert_eq!(
    ///        sus,
    ///        [Some(4),Some(3),Some(2),Some(4),Some(3),
    ///             Some(2),Some(1),Some(1)]
    ///    );
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
    /// this is raw intervals in Genome.seq
    pub fn find_intervals(&self, sus: Vec<Option<usize>>, kmer_size: usize) -> Vec<(usize, usize)> {
        let mut vec = Vec::new();
        let mut i: usize = 0;
        while i < self.length {
            let start = i; // the start of this interval
            let mut cur: usize = sus[i].unwrap_or(kmer_size + 1);
            // iterate until we do not find a unique substring
            while cur <= kmer_size // position is a SUN
                && (i + 1) < self.length // next position is not past the end of the genome
                && self.seq[i] != END_CHAR // current position is not the end of a contig
                && self.seq[i + 1] != END_CHAR
            // next pos is not the end of a contig
            {
                i += 1;
                cur = sus[i].unwrap_or(kmer_size + 1);
            }
            // add the interval to the return (+1 for bed interval)
            let end = i + 1;
            if end - start >= kmer_size {
                vec.push((start, end));
            }
            // increment to next start
            i += 1;
        }
        vec
    }

    /// Get back the chromosome possition from the flattened sequence
    fn convert_from_raw(
        &self,
        raw_intervals: Vec<(usize, usize)>,
    ) -> Vec<(&String, usize, usize, &[u8])> {
        let mut i = 0;
        let mut intervals = Vec::new();
        for (raw_start, raw_end) in raw_intervals {
            // increment our contig until we are in the correct interval
            while raw_start > self.ends[i] && raw_end > self.ends[i] {
                i += 1;
            }

            intervals.push((
                &self.names[i],
                raw_start - self.starts[i],
                raw_end - self.starts[i],
                &self.seq[raw_start..raw_end],
            ));
        }
        intervals
    }

    pub fn find_sun_intervals(&self, kmer_size: usize) -> Vec<(&String, usize, usize, &[u8])> {
        assert!(kmer_size > 1);
        let sus = Genome::get_shortest_subseq_size(&self.seq);
        let raw_intervals = self.find_intervals(sus, kmer_size);
        self.convert_from_raw(raw_intervals)
    }
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

pub fn validate_suns(
    genome: &Genome,
    intervals: &[(&String, usize, usize, &[u8])],
    kmer_size: usize,
) {
    let genome_str = std::str::from_utf8(&genome.seq).unwrap();
    println!("kmer size:{}\n\nSUN intervals:", kmer_size);
    let mut all_suns = Vec::new();
    // check that the SUNs only happen once and are valid
    for (chr, start, end, seq) in intervals {
        let sun_r = std::str::from_utf8(seq).unwrap();
        // print the result
        println!("{}\t{}\t{}\t{}", chr, start, end, sun_r);

        // break sun range into individual suns
        for i in 0..(sun_r.len() - kmer_size + 1) {
            let sun = &sun_r[i..i + kmer_size];
            all_suns.push(sun);
            // assert that the SUN only happens once
            assert_eq!(1, genome_str.matches(sun).count());
            // make sure no end characters made it into our SUN
            assert_eq!(0, sun.matches(END_CHAR_STR).count());
        }
    }
    // check that we found all sun kmers TODO
    for i in 0..(genome_str.len() - kmer_size) {
        let sun = &genome_str[i..i + kmer_size];
        // skip contig ends
        if sun.matches(END_CHAR_STR).count() > 0 {
            continue;
        }
        // break if we start looking at the reverse complement
        if i >= genome.length {
            break;
        }
        //println!("{}\tis this a sun?", sun);
        let count = genome_str.matches(sun).count();
        // assert that this kmer happens more than once or is a SUN
        assert!(count > 1 || all_suns.contains(&sun));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sun_finding() {
        let genome = Genome::from_file("test/test.fa");

        let kmer_size = 2;
        let intervals = genome.find_sun_intervals(kmer_size);
        validate_suns(&genome, &intervals, kmer_size);

        let kmer_size = 3;
        let intervals = genome.find_sun_intervals(kmer_size);
        validate_suns(&genome, &intervals, kmer_size);

        let kmer_size = 4;
        let intervals = genome.find_sun_intervals(kmer_size);
        validate_suns(&genome, &intervals, kmer_size);

        let kmer_size = 5;
        let intervals = genome.find_sun_intervals(kmer_size);
        validate_suns(&genome, &intervals, kmer_size);
    }
}
