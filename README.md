# rustybam

[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Test%20and%20Build/badge.svg)](https://github.com/mrvollger/rustybam/actions)
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Formatting/badge.svg)](https://github.com/mrvollger/rustybam/actions)
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Clippy/badge.svg)](https://github.com/mrvollger/rustybam/actions)
[![DOI](https://zenodo.org/badge/351639424.svg)](https://zenodo.org/badge/latestdoi/351639424)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/rustybam?color=green)](https://anaconda.org/bioconda/rustybam)
[![Downloads](https://img.shields.io/conda/dn/bioconda/rustybam?color=green)](https://anaconda.org/bioconda/rustybam)

<!--ts-->

- [rustybam](#rustybam)
  - [Usage](#usage)
    - [Available options and subcommands](#available-options-and-subcommands)
  - [Install](#install)
    - [conda](#conda)
    - [cargo](#cargo)
    - [Pre-complied binaries](#pre-complied-binaries)
    - [From source](#from-source)
  - [Examples](#examples)
    - [Manipulating PAFs and creating liftovers](#manipulating-pafs-and-creating-liftovers)
    - [Splitting up a fastx file](#splitting-up-a-fastx-file)
    - [Extract fasta from a genome with a bed file](#extract-fasta-from-a-genome-with-a-bed-file)
  - [TODO](#todo)

<!-- Added by: mrvollger, at: Tue Dec 21 21:29:21 PST 2021 -->

<!--te-->

## Usage

```
rustybam [OPTIONS] <SUBCOMMAND>
```

or

```
rb [OPTIONS] <SUBCOMMAND>
```

### Available options and subcommands

```
rustybam 0.1.23
Mitchell R. Vollger <mrvollger@gmail.com>
bioinformatics toolkit in rust

USAGE:
    rb [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -t, --threads <THREADS>    threads for decompression [default: 8]
    -v, --verbose              logging level
    -h, --help                 Print help information
    -V, --version              Print version information

SUBCOMMANDS:
    stats          Get percent identity stats from a sam/bam/cram or PAF
    bed-length     Count the number of bases in a bed file [aliases: bedlen, bl, bedlength]
    filter         Filter PAF records in various ways
    invert         Invert the target and query sequences in a PAF along with the CIGAR string
    liftover       Liftover target sequence coordinates onto query sequence using a PAF
    trim-paf       Trim paf records that overlap in query sequence [aliases: trim, tp]
    orient         Orient paf records so that most of the bases are in the forward direction
    break-paf      Break PAF records with large indels into multiple records (useful for
                   SafFire) [aliases: breakpaf, bp]
    paf-to-sam     Convert a PAF file into a SAM file. Warning, all alignments will be marked as
                   primary! [aliases: paftosam, p2s, paf2sam]
    fasta-split    Reads in a fasta from stdin and divides into files (can compress by adding
                   .gz) [aliases: fastasplit, fasplit]
    fastq-split    Reads in a fastq from stdin and divides into files (can compress by adding
                   .gz) [aliases: fastqsplit, fqsplit]
    get-fasta      Mimic bedtools getfasta but allow for bgzip in both bed and fasta inputs
                   [aliases: getfasta, gf]
    nucfreq        Get the frequencies of each bp at each position
    repeat         Report the longest exact repeat length at every position in a fasta
    suns           Extract the intervals in a genome (fasta) that are made up of SUNs
    help           Print this message or the help of the given subcommand(s)
```

## Install

### conda

```
mamba install -c bioconda rustybam
```

### cargo

```
cargo install rustybam
```

### Pre-complied binaries

Download from [releases](https://github.com/mrvollger/rustybam/releases) (may be slower than locally complied versions).

### From source

```

git clone https://github.com/mrvollger/rustybam.git
cd rustybam
cargo build --release

```

and the executable will be built here:

```

target/release/rustybam

```

There will also be an identical binary with the abbreviated name `rb`:

```

target/release/rb

```

## Examples

### Manipulating PAFs and creating liftovers

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

You can also use `rustybam breakpaf` to break up the paf records of indels above a certain size to
get more "miropeats" like intervals.

```

./rustybam breakpaf --max-size 1000 input.paf \
 | ./rustybam liftover \
 --bed <(printf "chr1\t0\t250000000\n") \
 | ./rustybam stats --paf \

> trimmed.stats.bed

```

> Yeah but how do I visualize the data?

Try out
[SafFire](https://mrvollger.github.io/SafFire/)!

### Splitting up a fastx file

Split a fasta file between `stdout` and two other files both compressed and uncompressed.

```shell
cat {input.fasta} | rustybam fasta-split two.fa.gz three.fa
```

Split a fastq file between `stdout` and two other files both compressed and uncompressed.

```shell
cat {input.fastq} | rustybam fastq-split two.fq.gz three.fq
```

### Extract fasta from a genome with a bed file

```
Mimic bedtools getfasta but allow for bgzip in both bed and fasta inputs

USAGE:
    rb get-fasta [OPTIONS] --bed <BED>

OPTIONS:
    -f, --fasta <FASTA>    fasta file to extract sequences from [default: -]
    -b, --bed <BED>        bed file of regions to extract
    -s, --strand           revcomp sequence if the strand is "-"
    -n, --name             add the name (4th column) to the header of the fasta output
    -h, --help             Print help information
    -V, --version          Print version information
```

## TODO

- [x] Finish implementing `trim-paf`.
- [x] Add a `bedtools getfasta` like operation that actually works with bgzipped input.
  - [ ] implement bed12/split
- [ ] Allow sam or paf for operations:
  - [x] make a sam header from a PAF file
  - [ ] convert sam record to paf record
  - [x] convert paf record to sam record
- [ ] Add `D4` for Nucfreq.
- [ ] Finish implementing `suns`.
