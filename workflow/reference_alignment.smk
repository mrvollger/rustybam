import os
import sys
import pandas as pd


df = pd.read_csv(config.get("tbl"), sep="\t")
df.asm = df.asm.map(os.path.abspath)
df["asm"] = df.asm.str.split(",")
df=df.explode("asm")
df['num'] = df.groupby(level=0).cumcount()+1
df.set_index( df["sample"] + "_" + df["num"].astype(str), inplace=True)

print(df)
wildcard_constraints:
    i="\d+"

def get_asm(wc):
    return df.loc[str(wc.sm)].asm


rule clean_query:
    input:
        query=get_asm,
    output:
        fasta="reference_alignment/{sm}.fasta",
    threads: 1
    run:
        import pysam
        out = open(output.fasta, "w+")
        seen = set()
        for idx, fasta in enumerate(input.query):
            for rec in pysam.FastxFile(fasta, persist=False):
                name = rec.name
                if rec.name in seen:
                    name = rec.name + f"_{idx+1}"
                seen.add(name)
                out.write(f">{name} {rec.comment}\n{rec.sequence}\n")
        out.close()


rule unimap:
    input:
        ref = os.path.abspath(config.get("ref")),
        query= get_asm,
    output:
        aln=pipe("reference_alignment/{sm}.sam"),
    log:
        "log/unimap.{sm}.log",
    benchmark:
        "log/unimap.{sm}.benchmark.txt"
    threads: config.get("aln_threads", 40)
    shell:
        """
        unimap -K 8G -t {threads} \
            -r 200000 -ax asm20 \
            --secondary=no --eqx -s 25000 \
                    {input.ref} {input.query} > {output.aln} \
                        2> {log}
        """


rule compress_sam:
    input:
        aln=rules.unimap.output.aln,
    output:
        aln="reference_alignment/bam/{sm}.bam",
        index="reference_alignment/bam/{sm}.bam.csi",
    threads: 8
    shell:
        """
        samtools view -@8 -u {input.aln} \
            | samtools sort -m 4G -@8 --write-index \
                 -o {output.aln}
        """


rule sam_to_paf:
    input:
        aln=rules.compress_sam.output.aln,
    output:
        paf="reference_alignment/paf/{sm}.paf.gz",
    threads: 8
    shell:
        "samtools view -@8 -h {input.aln} | paftools.js sam2paf - | gzip > {output.paf}"


rule paf_to_bed:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        bed="reference_alignment/bed/{sm}.bed",
    threads: 8
    shell:
        """
        printf "#t_nm\tt_st\tt_en\tq_nm\tq_st\tq_en\tstrand\tmapq\n" > {output.bed}
        gunzip -c {input.paf} \
            | awk -v OFS=$'\t' '{{print $6,$8,$9,$1,$3,$4,$5,$12}}' \
            >> {output.bed}
        """



rule reference_alignment:
    input:
        expand(rules.ra_sam_to_paf.output, sm=df.index),
        expand(rules.ra_paf_to_bed.output, sm=df.index),
    output:
    message: "Reference alignments complete"
