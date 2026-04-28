configfile: "config/config.yaml"

import pandas as pd
import glob

samples = pd.concat([
    pd.read_csv(sample_file, usecols=["sample", "library_ID", "short_R1_filename", "short_R2_filename", "long_filename"])
    for sample_file in glob.glob("samples/hybrid_assembly/*.csv")
])

rule all:
    input:
        expand("results/{library_ID}/{sample}/qc/fastp/{sample}_R1_clean.fastq.gz", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/qc/fastp/{sample}_R2_clean.fastq.gz", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/logs/{sample}_fastp_report.html", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/logs/{sample}_fastp_report.json", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/qc/Filtlong/{sample}_long_clean.fastq.gz", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/assembly/Autocycler/autocycler_out/consensus_assembly.fasta", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/qc/Polypolish/{sample}_polypolish.fasta", zip,
               library_ID=samples["library_ID"], sample=samples["sample"]),
        expand("results/{library_ID}/{sample}/stats/autocycler_assembly_metrics.tsv", zip,
               library_ID=samples["library_ID"], sample=samples["sample"])

rule fastp_qc_short_reads:
    threads: config["threads"]
    input:
        r1_reads=lambda wildcards: "raw_sequencing_reads/" + samples.loc[samples["sample"] == wildcards.sample, "short_R1_filename"].values[0],
        r2_reads=lambda wildcards: "raw_sequencing_reads/" + samples.loc[samples["sample"] == wildcards.sample, "short_R2_filename"].values[0]
    output:
        r1_clean="results/{library_ID}/{sample}/qc/fastp/{sample}_R1_clean.fastq.gz",
        r2_clean="results/{library_ID}/{sample}/qc/fastp/{sample}_R2_clean.fastq.gz",
        r1_unpaired="results/{library_ID}/{sample}/qc/fastp/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired="results/{library_ID}/{sample}/qc/fastp/{sample}_R2_unpaired.fastq.gz",
        html="results/{library_ID}/{sample}/logs/{sample}_fastp_report.html",
        json="results/{library_ID}/{sample}/logs/{sample}_fastp_report.json"
    log: "results/{library_ID}/{sample}/logs/{sample}_fastp.log"
    shell:
        """
        fastp \
            --in1 {input.r1_reads} --in2 {input.r2_reads} \
            --out1 {output.r1_clean} --out2 {output.r2_clean} \
            --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
              2> {log}
        """


rule Filtlong_qc_long_reads:
    threads: config["threads"]
    input:
        long_reads=lambda wildcards: "raw_sequencing_reads/" + samples.loc[samples["sample"] == wildcards.sample, "long_filename"].values[0]
    output:
        long_clean="results/{library_ID}/{sample}/qc/Filtlong/{sample}_long_clean.fastq.gz"
    log: "results/{library_ID}/{sample}/logs/{sample}_Filtlong.log"
    shell:
        """
        (filtlong --min_length 1000 --keep_percent 95 {input.long_reads} | gzip > {output.long_clean}) 2> {log}
        """

rule Autocycler_hybrid_assembly:
    threads: config["threads"]
    input:
        long_clean="results/{library_ID}/{sample}/qc/Filtlong/{sample}_long_clean.fastq.gz"
    output:
        assembly="results/{library_ID}/{sample}/assembly/Autocycler/autocycler_out/consensus_assembly.fasta"
    log: "results/{library_ID}/{sample}/logs/{sample}_Autocycler.log"
    shell:
        """
        #!/usr/bin/env bash

        # This script is a wrapper for running a fully-automated Autocycler assembly.

        # Usage:
        #   autocycler_full.sh <read_fastq> <threads> <jobs>

        # Copyright 2025 Ryan Wick (rrwick@gmail.com)
        # Licensed under the GNU General Public License v3.
        # See https://www.gnu.org/licenses/gpl-3.0.html.

        # Ensure script exits on error.
        set -e

        # Get arguments.
        reads={input.long_clean}                 # input reads FASTQ
        threads={threads}             # threads per job
        jobs=1                 # number of simultaneous jobs
        read_type="ont_r10"  # read type (default = ont_r10)

        # Input assembly jobs that exceed this time limit will be killed
        max_time="30h"

        # Validate input parameters
        if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
            echo "Usage: $0 <read_fastq> <threads> <jobs> [read_type]" 1>&2
            exit 1
        fi
        if [[ ! -f "$reads" ]]; then
            echo "Error: Input file '$reads' does not exist." 1>&2
            exit 1
        fi
        if (( threads > 128 )); then threads=128; fi  # Flye won't work with more than 128 threads

        case $read_type in
            ont_r9|ont_r10|pacbio_clr|pacbio_hifi) ;;
            *) echo "Error: read_type must be ont_r9, ont_r10, pacbio_clr or pacbio_hifi" 1>&2; exit 1 ;;
        esac

        genome_size=$(autocycler helper genome_size --reads "$reads" --threads "$threads")

        # Step 1: subsample the long-read set into multiple files
        mkdir -p results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler
        autocycler subsample --reads "$reads" --out_dir results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/subsampled_reads --genome_size "$genome_size" 2>> {log}

        # Step 2: assemble each subsampled file (myloasm and nextdenovo requre linux)
        mkdir -p results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/input_assemblies
        for assembler in raven myloasm miniasm flye metamdbg necat nextdenovo plassembler canu; do
            for i in 01 02 03 04; do
                autocycler helper "$assembler" --reads results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/subsampled_reads/sample_"$i".fastq --out_prefix results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/input_assemblies/"$assembler"_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
            done
        done

        # Give circular contigs from Plassembler extra clustering weight
        shopt -s nullglob
        for f in results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/input_assemblies/plassembler*.fasta; do
            perl -i -pe 's/circular=True/circular=True Autocycler_cluster_weight=3/' "$f"
        done

        # Give contigs from Canu and Flye extra consensus weight
        for f in results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/input_assemblies/canu*.fasta results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/input_assemblies/flye*.fasta; do
            perl -i -pe 's/^(>.*)$/$1 Autocycler_consensus_weight=2/' "$f"
        done
        shopt -u nullglob

        # Remove the subsampled reads to save space
        #rm results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/subsampled_reads/*.fastq

        # Step 3: compress the input assemblies into a unitig graph
        autocycler compress -i results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/input_assemblies -a results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/autocycler_out 2>> {log}

        # Step 4: cluster the input contigs into putative genomic sequences
        autocycler cluster -a results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/autocycler_out 2>> {log}

        # Steps 5 and 6: trim and resolve each QC-pass cluster
        for c in results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/autocycler_out/clustering/qc_pass/cluster_*; do
            autocycler trim -c "$c" 2>> {log}
            autocycler resolve -c "$c" 2>> {log}
        done

        # Step 7: combine resolved clusters into a final assembly
        autocycler combine -a results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/autocycler_out -i results/{wildcards.library_ID}/{wildcards.sample}/assembly/Autocycler/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> {log}
        """

rule polypolish_short_read_polishing:
    threads: config["threads"]
    input:
        assembly="results/{library_ID}/{sample}/assembly/Autocycler/autocycler_out/consensus_assembly.fasta",
        r1_clean="results/{library_ID}/{sample}/qc/fastp/{sample}_R1_clean.fastq.gz",
        r2_clean="results/{library_ID}/{sample}/qc/fastp/{sample}_R2_clean.fastq.gz",
    output:
        r1_alignment="results/{library_ID}/{sample}/qc/Polypolish/{sample}_R1_alignment.sam",
        r2_alignment="results/{library_ID}/{sample}/qc/Polypolish/{sample}_R2_alignment.sam",
        r1_alignment_filtered="results/{library_ID}/{sample}/qc/Polypolish/{sample}_R1_alignment_filtered.sam",
        r2_alignment_filtered="results/{library_ID}/{sample}/qc/Polypolish/{sample}_R2_alignment_filtered.sam",
        polished_consensus_assembly="results/{library_ID}/{sample}/qc/Polypolish/{sample}_polypolish.fasta"
    log: "results/{library_ID}/{sample}/logs/{sample}_Polypolish.log"
    shell:
        """
        bwa-mem2 index {input.assembly}
        bwa-mem2 mem -t {threads} -a {input.assembly} {input.r1_clean} > {output.r1_alignment}
        bwa-mem2 mem -t {threads} -a {input.assembly} {input.r2_clean} > {output.r2_alignment}
        polypolish filter --in1 {output.r1_alignment} --in2 {output.r2_alignment} --out1 {output.r1_alignment_filtered} --out2 {output.r2_alignment_filtered}
        polypolish polish {input.assembly} {output.r1_alignment_filtered} {output.r2_alignment_filtered} > {output.polished_consensus_assembly}
        """

rule assembly_stats:
    threads: config["threads"]
    input:
        assembly_dir="results/{library_ID}/{sample}/assembly/Autocycler",
        consensus="results/{library_ID}/{sample}/assembly/Autocycler/autocycler_out/consensus_assembly.fasta",
        polished_consensus_assembly="results/{library_ID}/{sample}/qc/Polypolish/{sample}_polypolish.fasta"
    output:
        assembly_stats="results/{library_ID}/{sample}/stats/autocycler_assembly_metrics.tsv"
    shell:
        """
        assembly_date=$(date +%Y%m%d)
        assembler_version=$(autocycler --version)

        echo -e "path\tinput_read_count\tinput_read_bases\tinput_read_n50\tpass_cluster_count\tfail_cluster_count\toverall_clustering_score\tuntrimmed_cluster_size\tuntrimmed_cluster_distance\ttrimmed_cluster_size\ttrimmed_cluster_median\ttrimmed_cluster_mad\tconsensus_assembly_bases\tconsensus_assembly_unitigs\tconsensus_assembly_fully_resolved\tlibrary_ID\tsample\tassembly_date\tassembler" > {output.assembly_stats}
        autocycler table -a {input.assembly_dir} -n {input.polished_consensus_assembly} | \
            awk -v library_ID={wildcards.library_ID} -v sample={wildcards.sample} -v assembly_date="$assembly_date" -v assembler="$assembler_version" \
                '{{print $0"\t"library_ID"\t"sample"\t"assembly_date"\t"assembler}}' >> {output.assembly_stats}
        """