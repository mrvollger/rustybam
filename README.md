# rustybam
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Test%20and%20Build/badge.svg)](https://github.com/mrvollger/rustybam/actions) 
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Formatting/badge.svg)](https://github.com/mrvollger/rustybam/actions) 
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Clippy/badge.svg)](https://github.com/mrvollger/rustybam/actions) 

# Install
It is easy just make sure you have rust installed and then:
```
git clone https://github.com/mrvollger/rustybam.git
cd rustybam 
cargo build --release 
```
and the executable will be built here:
```
target/release/rustybam 
```
# Examples 
> I have a `PAF` and I want to subset it for just a particular region in the reference. 

With `rustybam` its easy:
```
./rustybam liftover \
    --bed <(printf "chr1\t0\t250000000\n") \
    input.paf > trimmed.paf 
```
> But I also want the alignment statistics for the region. 

No problem, `rustybam liftover` does not just trim the coordinates but also the CIGAR
so it is ready for `rustybam stats`:
```
./rustybam liftover \
    --bed <(printf "chr1\t0\t250000000\n") \
    input.paf \
    | ./rustybam stats --paf \
    > trimmed.stats.bed 
```
> Okay, but Evan asked for an "align slider" so I need to realign in chunks. 

No need, just make your `bed` query to `rustybam liftoff` a set of sliding windows 
and it will do the rest. 
```
./rustybam liftover \
    --bed <(bedtools makewindows -w 100000 \
                <(printf "chr1\t0\t250000000\n")
            ) \
    input.paf \
    | ./rustybam stats --paf \
    > trimmed.stats.bed 
```
> Yeah but how do I visualize the data?

Try out 
[Miropeats in D3](https://mrvollger.github.io/D3Miropeats/)!

# Usage
```
./rustybam 0.0.1
Mitchell R. Vollger's alignment utilities

USAGE:
    rustybam [OPTIONS] <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -t, --threads <threads>    Number of threads to use for decompressing

SUBCOMMANDS:
    help        Prints this message or the help of the given subcommand(s)
    liftover    liftover target sequence coordinates onto query sequence using a PAF
    nucfreq     Get the freqs of each bp at each position.
    stats       Get percent identity stats from a sam/bam/cram or PAF (add --paf)
    suns        Extract the intervals in a genome (fasta) that are made up of SUNs
```
