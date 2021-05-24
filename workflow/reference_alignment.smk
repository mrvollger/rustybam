import os
import sys
import pandas as pd


df = pd.read_csv(config.get("tbl"), sep="\t")
df.asm = df.asm.map(os.path.abspath)
df["asm"] = df.asm.str.split(",")
df = df.explode("asm")
df["num"] = df.groupby(level=0).cumcount() + 1
df.set_index(df["sample"] + "_" + df["num"].astype(str), inplace=True)


wildcard_constraints:
    i="\d+",


def get_asm(wc):
    return df.loc[str(wc.sm)].asm


def get_ref(wc):
    return config.get("ref")[wc.ref]


def get_fai(wc):
    return config.get("ref")[wc.ref] + ".fai"


rule unimap_index:
    input:
        ref=get_ref,
    output:
        umi="reference_alignment/{ref}/{ref}.umi",
    threads: 8
    conda:
        "envs/environment.yml"
    shell:
        "unimap -t {threads} -ax asm20 -d {output.umi} {input.ref}"


rule unimap:
    input:
        ref=rules.unimap_index.output.umi,
        query=get_asm,
    output:
        aln=pipe("reference_alignment/{ref}/{sm}.sam"),
    log:
        "log/unimap.{ref}_{sm}.log",
    benchmark:
        "log/unimap.{ref}_{sm}.benchmark.txt"
    conda:
        "envs/environment.yml"
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
        aln="reference_alignment/{ref}/bam/{sm}.bam",
        index="reference_alignment/{ref}/bam/{sm}.bam.csi",
    threads: 1 # dont increase this, it will break things randomly 
    conda:
        "envs/environment.yml"
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
        paf="reference_alignment/{ref}/paf/{sm}.paf",
    conda:
        "envs/environment.yml"
    shell:
        "samtools view -h {input.aln} | paftools.js sam2paf - > {output.paf}"


rule paf_to_bed:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        bed="reference_alignment/{ref}/bed/{sm}.bed",
    threads: 8
    conda:
        "envs/environment.yml"
    params:
        rb=config["rb"],
    shell:
        """
        {params.rb} stats --paf {input.paf} > {output.bed}
        """


rule bed_to_pdf:
    input:
        bed="reference_alignment/{ref}/bed/{sm}_1.bed",
        bed2="reference_alignment/{ref}/bed/{sm}_2.bed",
    output:
        pdf="reference_alignment/{ref}/pdf/ideogram.{sm}.pdf",
    threads: 1
    conda:
        "envs/environment.yml"
    params:
        smkdir=config["smkdir"],
    shell:
        """
        Rscript {params.smkdir}/scripts/ideogram.R \
          --asm {input.bed} \
          --asm2 {input.bed2} \
          --plot {output.pdf}
        """


rule query_ends:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        bed=temp("reference_alignment/{ref}/ends/tmp.{sm}.bed"),
    params:
        smkdir=config["smkdir"],
    conda:
        "envs/environment.yml"
    threads: 1
    shell:
        """
        {params.smkdir}/scripts/ends_from_paf.py \
          --minwidth 50000 \
          --width 1000 \
          {input.paf} > {output.bed}
        """


rule find_contig_ends:
    input:
        paf=rules.sam_to_paf.output.paf,
        bed=rules.query_ends.output.bed,
    output:
        bed="reference_alignment/{ref}/ends/{sm}.bed",
    threads: 1
    conda:
        "envs/environment.yml"
    params:
        rb=config["rb"],
    shell:
        """
        {params.rb} liftover --largest --qbed --bed {input.bed} {input.paf} \
          | {params.rb} stats --paf --qbed \
          > {output.bed}
        """


rule collect_contig_ends:
    input:
        beds=expand(rules.find_contig_ends.output.bed, sm=df.index, ref="{ref}"),
    output:
        bed="reference_alignment/{ref}/ends/all.ends.bed",
    threads: 1
    conda:
        "envs/environment.yml"
    shell:
        """
        head -n 1 {input.beds[0]} > {output.bed}
        cat {input.beds} \
          | grep -v "^#" \
          | bedtools sort -i - \
          >> {output.bed}
        """


rule windowed_ends:
    input:
        fai=get_fai,
        bed=rules.collect_contig_ends.output.bed,
    output:
        bed="reference_alignment/{ref}/ends/windowed.all.ends.bed",
    threads: 1
    conda:
        "envs/environment.yml"
    shell:
        """
        bedtools intersect -wa -wb -header \
          -a <(printf "#chr\tstart\tend\n" ; bedtools makewindows -w 1000000 -g {input.fai} ) \
          -b {input.bed} \
          > {output.bed}
        header=$(head -n1 {input.bed})
        sed -i " 1 s/$/\t$header/" {output.bed}
        """


rule pre_end_content:
    input:
        ref=get_ref,
        fai=get_fai,
    output:
        allbed="reference_alignment/{ref}/ends/all.nuc.content.bed",
    threads: 1
    conda:
        "envs/environment.yml"
    shell:
        """
        bedtools nuc \
          -fi {input.ref} \
          -bed <(bedtools makewindows -s 100 -w 1000 -g {input.fai} ) \
          > {output.allbed}
        """


rule end_content:
    input:
        bed=rules.collect_contig_ends.output.bed,
        allbed=rules.pre_end_content.output.allbed,
        fai=get_fai,
    output:
        bed="reference_alignment/{ref}/ends/all.ends.nuc.content.bed",
    threads: 1
    conda:
        "envs/environment.yml"
    shell:
        """
        bedtools intersect -header -u -a {input.allbed} \
          -b <(bedtools slop -b 10000 -g {input.fai} -i {input.bed}) \
          > {output.bed}
        """


rule reference_alignment:
    input:
        expand(rules.collect_contig_ends.output, ref=config.get("ref").keys()),
        expand(rules.end_content.output, ref=config.get("ref").keys()),
        expand(rules.windowed_ends.output, ref=config.get("ref").keys()),
        expand(rules.ra_sam_to_paf.output, sm=df.index, ref=config.get("ref").keys()),
        expand(
            rules.bed_to_pdf.output,
            sm=df["sample"].str.strip(),
            ref=config.get("ref").keys(),
        ),
        expand(rules.ra_paf_to_bed.output, sm=df.index, ref=config.get("ref").keys()),
        expand(rules.find_contig_ends.output, sm=df.index, ref=config.get("ref").keys()),
    message:
        "Reference alignments complete"
