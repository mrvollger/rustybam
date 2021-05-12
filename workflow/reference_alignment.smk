import os
import sys
import pandas as pd


df = pd.read_csv(config.get("tbl"), sep="\t")
df.asm = df.asm.map(os.path.abspath)
df["asm"] = df.asm.str.split(",")
df = df.explode("asm")
df["num"] = df.groupby(level=0).cumcount() + 1
df.set_index(df["sample"] + "_" + df["num"].astype(str), inplace=True)
shell("which python")

wildcard_constraints:
    i="\d+",

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


rule unimap_index:
    input:
        ref=os.path.abspath(config.get("ref")),
    output:
        umi="reference_alignment/ref.umi",
    threads: 8
    shell:
        "unimap -t {threads} -ax asm20 -d {output.umi} {input.ref}"


rule unimap:
    input:
        ref=rules.unimap_index.output.umi,
        query=get_asm,
    output:
        aln=pipe("reference_alignment/{sm}.sam"),
    log:
        "log/unimap.{sm}.log",
    benchmark:
        "log/unimap.{sm}.benchmark.txt"
    threads: config.get("aln_threads", 4)
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
    threads: 1  # dont increase this, it will break things randomly 
    shell:
        """
        samtools view -u {input.aln} \
            | samtools sort -m 20G --write-index \
                 -o {output.aln}
        """


rule sam_to_paf:
    input:
        aln=rules.compress_sam.output.aln,
    output:
        paf="reference_alignment/paf/{sm}.paf",
    shell:
        "samtools view -h {input.aln} | paftools.js sam2paf - > {output.paf}"


rule paf_to_bed:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        bed="reference_alignment/bed/{sm}.bed",
    threads: 8
    params:
      rb = config["rb"]
    shell:
        """
        {params.rb} stats --paf {input.paf} > {output.bed}
        """


rule query_ends:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        bed=temp("reference_alignment/ends/tmp.{sm}.bed"),
    params:
      smkdir = config["smkdir"]
    threads: 1
    shell: "which python; {params.smkdir}/scripts/ends_from_paf.py --width 50000 {input.paf} > {output.bed}"


rule find_contig_ends:
    input:
        paf=rules.sam_to_paf.output.paf,
        bed=rules.query_ends.output.bed,
    output:
        bed="reference_alignment/ends/{sm}.bed",
    threads: 1
    params:
      rb = config["rb"]
    shell:
        """
        {params.rb} liftover --qbed --bed {input.bed} {input.paf} \
          | {params.rb} stats --paf \
          > {output.bed}
        """


rule collect_contig_ends:
    input:
        beds=expand(rules.find_contig_ends.output.bed, sm=df.index),
    output:
        bed="reference_alignment/ends/all.ends.bed",
    threads: 1
    shell:
        """
        head -n 1 {input.beds[0]} > {output.bed}
        cat {input.beds} | grep -v "^#" >> {output.bed}
        """


rule reference_alignment:
    input:
        rules.collect_contig_ends.output,
        expand(rules.ra_sam_to_paf.output, sm=df.index),
        expand(rules.ra_paf_to_bed.output, sm=df.index),
        expand(rules.find_contig_ends.output, sm=df.index),
    message:
        "Reference alignments complete"
