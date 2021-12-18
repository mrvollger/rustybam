# rustybam

[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Test%20and%20Build/badge.svg)](https://github.com/mrvollger/rustybam/actions)
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Formatting/badge.svg)](https://github.com/mrvollger/rustybam/actions)
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Clippy/badge.svg)](https://github.com/mrvollger/rustybam/actions)
[![DOI](https://zenodo.org/badge/351639424.svg)](https://zenodo.org/badge/latestdoi/351639424)
[![](https://anaconda.org/bioconda/rustybam/badges/downloads.svg)](https://anaconda.org/bioconda/rustybam)
[![](https://anaconda.org/bioconda/rustybam/badges/installer/conda.svg)](https://anaconda.org/bioconda/rustybam)

## Install

It is easy just make sure you have rust and cargo installed and then:

```
cargo install rustybam
```

Alternatively to build from the latest source code:

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

## rustybam getfasta

```
rb-getfasta

Mimic bedtools getfasta but allow for bgzip in both bed and fasta inputs.

USAGE:
    rb getfasta [FLAGS] [OPTIONS] --bed <bed> --fi <fasta>

FLAGS:
    -h, --help       Print help information
    -n, --name       Add the name (4th column) to the header of the fasta output
    -s, --strand     revcomp sequence if the strand is reverse (-)
    -V, --version    Print version information

OPTIONS:
    -b, --bed <bed>            Bed file of regions to extract
    -f, --fi <fasta>           Fasta file to extract sequences from
    -t, --threads <threads>    Number of threads to use for decompressing
```

### More details on `liftover`

This is a function for lifting over coordinates from a reference (`--bed`) to a query using a `PAF` file from `minimap2` or `unimap`. `minimap2` (or `unimap`) must be run with `-c --eqx` and the output format must be `PAF` or else the liftover is not possible. Note, you can use `paftools.js sam2paf` to convert `sam/bam/cram` data to `paf` format for `liftover`.

The returned file is a `PAF` file that is trimmed to the regions in the bed file. Even the cigar in the returned PAF file is trimmed so it can be used downstream! Additionally, a tag with the format `id:Z:<>` is added to the `PAF` where `<>` is either the 4th column of the input bed file or if not present `chr_start_end`.

Want to liftover from the query to the reference? No problem, just pass the `-q` flag. Note, that this will make the query in the input `PAF` the target in the output `PAF`.

```
rustybam-liftover
liftover target sequence coordinates onto query sequence using a PAF

USAGE:
    rustybam liftover [FLAGS] [OPTIONS] --bed <bed> [paf]

ARGS:
    <paf>    PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero
             unless the cigar uses =X.

FLAGS:
    -h, --help       Prints help information
    -q, --qbed       The bed contains query coordinates to be lifted (note the query in the original
                     PAF will become the target in the output)
    -V, --version    Prints version information

OPTIONS:
    -b, --bed <bed>            Bed file of regions to liftover
    -t, --threads <threads>    Number of threads to use for decompressing
```

### Alignment workflow

This repository also includes a snakemake workflow for aligning genome assemblies to a reference. The config file `config/config.yaml` should be configured to your use case and the file `config/table.asm.tbl` should be configured for your input assemblies. The `Snakefile` can be found under the `workflow` directory.

Once you have modified `table.asm.tbl` to have your assemblies and `config.yaml` with your references you can run the snakemake with:

```
snakemake all --cores 160 --use-conda  {any extra snakemake options}
```

## TODO

- Add a `bedtools getfasta` like operation that actually works with bgzipped input.
  - Done, but no bed12/split implementation yet.
