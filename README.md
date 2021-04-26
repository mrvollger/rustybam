# rustybam
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Test%20and%20Build/badge.svg)](https://github.com/mrvollger/rustybam/actions) 
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Formatting/badge.svg)](https://github.com/mrvollger/rustybam/actions) 
[![Actions Status](https://github.com/mrvollger/rustybam/workflows/Clippy/badge.svg)](https://github.com/mrvollger/rustybam/actions) 

# Usage 
```
./rustybam 
Mitchell R. Vollger's bam utilities.

USAGE:
    rustybam [OPTIONS] <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -t, --threads <threads>    Number of threads to use for uncompressing

SUBCOMMANDS:
    help       Prints this message or the help of the given subcommand(s)
    nucfreq    Get the freqs of each bp at each postion.
    stats      Get precent identity stats from a sam/bam/cram (PAF coming soon).
```

## NucFreq
```
./rustybam nucfreq --help
Get the freqs of each bp at each postion.

USAGE:
    rustybam nucfreq [FLAGS] [OPTIONS] [BAM]

ARGS:
    <BAM>    sam/bam/cram file

FLAGS:
    -h, --help       Prints help information
    -s, --small      smaller output format.
    -V, --version    Prints version information

OPTIONS:
    -b, --bed <bed>            Print nucfreq info from regions in the bed file, output is optionally
                               tagged using the 4th column.
    -r, --region <region>      Print nucfreq info from the input region e.g "chr1:1-1000".
    -t, --threads <threads>    Number of threads to use for uncompressing

```

## Stats (samIdentity)
```
./rustybam stats --help
Get precent identity stats from a sam/bam/cram (PAF coming soon).

USAGE:
    rustybam stats [FLAGS] [OPTIONS] [BAM]

ARGS:
    <BAM>    sam/bam/cram file

FLAGS:
    -h, --help       Prints help information
    -q, --qbed       Print query coordinates first
    -V, --version    Prints version information

OPTIONS:
    -t, --threads <threads>    Number of threads to use for uncompressing
```

## Liftover 
This is a function for lifting over coordinates from a reference (`-bed`) 
  to a query using a PAF file from `minimap2` or `unimap`. 

The return file is a PAF file that is trimmed to the regions in the bed file. 
Even the cigar in the returned PAF file is trimmed so it can be used downstream!

Want to liftover from the query to the reference? 
  No problem, just pass the `-q` flag. Note, that this will make
  the query in the input PAF the target in the output PAF. 

```
rustybam-liftover 
liftover target sequence coordinates onto query sequence using a paf

USAGE:
    rustybam liftover [FLAGS] [OPTIONS] --bed <bed> [paf]

ARGS:
    <paf>    PAF file from minimap2 or unimap. Must have the cg tag

FLAGS:
    -h, --help        Prints help information
    -q, --querybed    The bed contains query coordiantes to be lifted (note the query in the
                      original PAF will become the target in the output)
    -V, --version     Prints version information

OPTIONS:
    -b, --bed <bed>            Bed file of regions to liftover
    -t, --threads <threads>    Number of threads to use for uncompressing
```

