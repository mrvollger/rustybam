```bash
rb trim-paf .test/asm_small.paf `#trims back alignments that align the same query sequence more than once` \
    | rb break-paf --max-size 100 `#breaks the alignment into smaller pieces on indels of 100 bases or more` \
    | rb orient `#orients each contig so that the majority of bases are forward aligned` \
    | rb liftover --bed <(printf "chr22\t12000000\t13000000\n") `#subsets and trims the alignment to 1 Mbp of chr22.` \
    | rb filter --paired-len 10000 `#filters for query sequences that have at least 10,000 bases aligned to a target across all alignments.` \
    | rb stats --paf `#calculates statistics from the trimmed paf file` \
    | less -S
```
