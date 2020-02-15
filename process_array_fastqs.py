#!/usr/bin/env python

"""
Trim and align fastqs from an array library to a bowtie2 reference.
Does NOT remove PCR duplicates. 
Creates a stranded bed file with original cluster_ID in 'name' column

Note: Python 3

Ben Ober-Reynolds
Stanford University
"""

import pandas as pd
import os
import sys
import glob
#import json


include: "process_array_fastqs_config.py" 


# Make metadata file if it doesn't exist
if not os.path.exists(METADATA_FILE): make_meta(METADATA_FILE)

# Add execution directory to path if it isn't already there
if EXE_DIR not in sys.path: os.environ["PATH"] = EXE_DIR + os.pathsep + os.environ["PATH"]

# Load metadata
metadata = pd.read_table(METADATA_FILE, index_col = False)

sample_labels = metadata.Name.tolist()


rule all:
    input:
        # Per sample output files
        # These are listed in the order generated
        expand("output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html", sample_label = sample_labels),
        expand("output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html", sample_label = sample_labels),
        expand("output/beds/{sample_label}.cluster.bed.gz", sample_label = sample_labels),

        # Pooled stats, etc.
        "output/bams/qc/compiled_flagstats.txt",
        "output/plots/qc/insert_size/all_sample_fragment_dist.pdf",
        # Pooled plots:
        "output/plots/qc/compiled_flagstats.pdf",

    output:
        "snakeATAC.txt"
    shell:
        "echo $(date) > {output};"
        "echo snake make stuff"


"""
Trim Nextera adapters using Skewer 
Version 0.2.2
"""

rule trim_adapters_skewer:
    input:
        left = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read1"].values[0]),
        right = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read2"].values[0])
    output:
        temp_left_cat = temp("output/fastqs/{sample_label}_skewer_R1.fastq.gz"),
        temp_right_cat = temp("output/fastqs/{sample_label}_skewer_R2.fastq.gz"),
        left = temp("output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz"),
        right = temp("output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"),
        log = "output/fastqs/qc/{sample_label}.skewer.log"
    params:
        # Cluster params
        rule_name = "trim_adapters_skewer",
        run_time = "2:30:00",
        cores = "4",
        memory = "6GB",
        job_name = "trimming_{sample_label}"
    benchmark:
        "benchmarks/trimming/{sample_label}.txt"
    threads: 4
    shell:
        "cat {input.left} > {output.temp_left_cat};" # if there are multiple files to be combined
        "cat {input.right} > {output.temp_right_cat};"
        "skewer \
        -x CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
        -y CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -m pe -t {threads} -f sanger\
        {output.temp_left_cat} {output.temp_right_cat} \
        -o output/fastqs/trimmed/{wildcards.sample_label} -z --quiet;"
        "mv output/fastqs/trimmed/{wildcards.sample_label}-trimmed-pair1.fastq.gz {output.left};"
        "mv output/fastqs/trimmed/{wildcards.sample_label}-trimmed-pair2.fastq.gz {output.right};"
        "mv output/fastqs/trimmed/{wildcards.sample_label}-trimmed.log {output.log};"


"""
Run fastQC on the trimmed and untrimmed fastqs to get some information about
potential problems with fastqs
"""
rule fastqc_unmapped_trimmed:
    input:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right = "output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"
    output:
        "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html",
        "output/fastqs/qc/{sample_label}_R2_trimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.zip",
        "output/fastqs/qc/{sample_label}_R2_trimmed_fastqc.zip"
    params:
        # Cluster params
        rule_name = "fastqc_unmapped_trimmed",
        run_time="00:59:00",
        cores="1",
        memory="6GB",
        job_name="fastqc_um_tr_{sample_label}"
    benchmark: 
        "benchmarks/fastqc/{sample_label}_trim.txt"
    shell:
        "fastqc {input.left} {input.right} --outdir=" + "output/fastqs/qc/"


rule fastqc_unmapped_untrimmed:
    input:
        left = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read1"].values[0]),
        right = lambda wildcards: glob.glob(metadata.loc[metadata.Name == wildcards.sample_label]["Read2"].values[0]),
    output:
        temp_left = temp("output/fastqs/qc/{sample_label}_R1_untrimmed.fastq.gz"),
        temp_right = temp("output/fastqs/qc/{sample_label}_R2_untrimmed.fastq.gz"),
        lh = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html",
        rh = "output/fastqs/qc/{sample_label}_R2_untrimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        lz = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.zip",
        rz = "output/fastqs/qc/{sample_label}_R2_untrimmed_fastqc.zip"
    params:
        # Cluster params
        rule_name = "fastqc_unmapped_untrimmed",
        run_time="00:59:00",
        cores="1",
        memory="6GB",
        job_name="fastqc_um_ut_{sample_label}"
    benchmark: 
        "benchmarks/fastqc/{sample_label}_untrim.txt"
    run:
        shell("cat {input.left} > {output.temp_left}"),
        shell("cat {input.right} > {output.temp_right}"),
        shell("fastqc {output.temp_left} {output.temp_right} --outdir=output/fastqs/qc/;")


"""
Map trimmed reads using Bowtie2
Version 2.2.6
Excludes mates separated by more than 2000 bp
Sorts and indexes the bam file afterwards using samtools
For info on sam (sequence alignment map) and bam (binary of sam):
https://training.h3abionet.org/postgraduate_workshop_2014/wp-content/uploads/2014/04/H3ABioNet_2014_NGS_8_SamFormat.pdf
"""
rule run_bowtie:
    input:
        # Adding the '.1.bt2' is necessary for snakemake to recognize the file
        idx = REFERENCE_FILE + ".1.bt2",
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        right ="output/fastqs/trimmed/{sample_label}_R2_trimmed.fastq.gz"
    output:
        bam = temp("output/bams/unprocessed/{sample_label}.bam"),
        idx = temp("output/bams/unprocessed/{sample_label}.bam.bai")
    params:
        # Cluster params
        rule_name = "run_bowtie",
        run_time = "8:00:00",
        cores = "8",
        memory = "20GB",
        job_name = "bwt2_{sample_label}"
    benchmark: "benchmarks/bowtie/{sample_label}.txt"
    threads: 8
    shell: 
        # -X 2000 # prevents mates separated by a lot 
        "bowtie2 \
        -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50\
        -X 2000 --threads {threads} \
        --rg-id {wildcards.sample_label} \
        --rg 'SM:{wildcards.sample_label}' \
        -x " + REFERENCE_FILE + " -1 {input.left} -2 {input.right} \
        | samtools view -b -S - \
        | samtools sort -o output/bams/unprocessed/{wildcards.sample_label}.bam -; " # -n flag makes pairs get sorted together
        "samtools index output/bams/unprocessed/{wildcards.sample_label}.bam; "



"""
flagstats calculated with SAMTools
"""
rule calc_flagstats:
    input:
        bam = "output/bams/unprocessed/{sample_label}.bam",
        idx = "output/bams/unprocessed/{sample_label}.bam.bai"
    output:
        "output/bams/qc/flagstats/{sample_label}.flagstat.txt" 
    params:
        # Cluster params
        rule_name = "calc_flagstats",
        run_time="00:10:00",
        cores="1",
        memory="3GB",
        job_name="flagstat_{sample_label}"
    shell:
        "samtools flagstat {input.bam} | awk '{{print \"{wildcards.sample_label}\\t\" $0}}' > {output};"


"""
idxstats calculated with SAMTools
"""
rule calc_idxstats:
    input:
        bam = "output/bams/unprocessed/{sample_label}.bam",
        idx = "output/bams/unprocessed/{sample_label}.bam.bai"
    output:
        "output/bams/qc/idxstats/{sample_label}.idxstats.txt" 
    params:
        # Cluster params
        rule_name = "calc_idxstats",
        run_time="00:10:00",
        cores="1",
        memory="1GB",
        job_name="idxstats_{sample_label}"
    shell:
        "samtools idxstats {input.bam} | awk '{{print \"{wildcards.sample_label}\\t\" $0}}' > {output};"


"""
compile SAMTools flagstats of all samples into one table.
"""
rule plot_flagstats:
    input:
        expand("output/bams/qc/flagstats/{sample_label}.flagstat.txt", sample_label=sample_labels)
    output:
        table = "output/bams/qc/compiled_flagstats.txt",
        pdf = "output/plots/qc/compiled_flagstats.pdf"
    params:
        # Cluster params
        rule_name = "plot_flagstats",
        run_time="00:10:00",
        cores="1",
        memory="1GB",
        job_name="plot_flagstat"
    shell:
        (
            "awk 'BEGIN {{OFS = \"\\t\"; print \"sample_label\",\"total\",\"secondary\","
            "\"supplementary\",\"duplicates\",\"mapped\",\"paired\",\"read1\",\"read2\","
            "\"proper_pair\",\"both_mapped\",\"singletons\",\"separate_chr\",\"separate_chr_mapq_above5\"}} "
            "FNR == 1 && NR != 1 {{print \"\"}} FNR == 1 {{printf $1}} {{printf \"\\t\" $2 }} "
            "END {{print \"\"}} ' {input} > {output.table};"
            "Rscript --vanilla {ARRAY_TOOLS}/qc_boxplot.R {output.table} read_count {output.pdf}"
        )


"""
Filter bams of low quality or unmapped reads
"""
rule filter_bams:
    input: 
        bam = rules.run_bowtie.output.bam,
        idx = rules.run_bowtie.output.idx
    output: 
        bam = "output/bams/filtered/{sample_label}.filtered.bam",
        idx = "output/bams/filtered/{sample_label}.filtered.bam.bai"
    params:
        # Cluster params
        rule_name = "filter_bams",
        run_time="05:00:00",
        cores="1",
        memory="8000",
        job_name="filter_bam_{sample_label}",
        # Rule params
        mapq_threshold="20"
    threads: 1
    run:
        # -F 1804: exclude flag, exludes unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
        # See for explaination: https://broadinstitute.github.io/picard/explain-flags.html
        # -f 2: flags to require, properly aligned
        # -q 30: exlude low MAPQ, set as parameter to adjust
        if BLACKLIST is None:
            shell(
                "samtools view -F 1804 -f 2 -q {params.mapq_threshold} -b {input.bam} > {output.bam}; \
                samtools index {output.bam}; "
                )
        else:
            shell(
                "samtools view -F 1804 -f 2 -q {params.mapq_threshold} -b {input.bam} \
                | bedtools intersect -v -abam - -b " + BLACKLIST + " -wa > {output.bam}; \
                samtools index {output.bam}; "
                )


"""
Plot the insert size distribution of each sample
"""
rule plot_insert_sizes:
    input:
        bams = expand("output/bams/filtered/{sample_label}.filtered.bam", sample_label=sample_labels)
    output:
        all_samp_ins = "output/plots/qc/insert_size/all_sample_fragment_dist.pdf",
    params:
        # Cluster params
        rule_name = "plot_insert_sizes",
        run_time="02:00:00",
        cores="1",
        memory="50GB",
        job_name="plot_insert_size",
        # Rule params
        bam_dir="output/bams/filtered",
        plot_dir="output/plots/qc/insert_size"
    threads: len(sample_labels)
    shell:
        (
            "Rscript --vanilla {ARRAY_TOOLS}/plot_fragment_length_dists.R {params.bam_dir} {params.plot_dir} {threads}"
        )


"""
Create an insertion bed file for each sample. 
Also adjusts the insertion site based on how the transposase sits on the DNA
"""
rule make_cluster_bed:
    input:
        bam = rules.filter_bams.output.bam,
        idx = rules.filter_bams.output.idx
    output:
        bed = "output/beds/{sample_label}.cluster.bed.gz",
        fasta = "output/beds/{sample_label}.cluster.fasta"
    params:
        # Cluster params
        rule_name = "make_insertion_bed",
        run_time="2:00:00",
        cores="1",
        memory="10GB",
        job_name="bam2bed_{sample_label}"
    threads: 1
    benchmark: "benchmarks/make_bed/{sample_label}_bam2bed.txt"
    shell:
        # Create bed file and then attach cluster information
        # Additionally, retrieve fastas for the actual array fragments
        (
            
            "Rscript --vanilla {ARRAY_TOOLS}/make_fragment_bed.R {input.bam} {output.bed}; "
            "bedtools getfasta -fi {REFERENCE_FILE}.fa -bed {output.bed} -fo {output.fasta} -name -s" 
        )  
  



