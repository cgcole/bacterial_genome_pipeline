#!/usr/bin/env bash

usage() {
    echo "Usage: $0 --threads <int> --samples <file>"
    echo "  --threads    Number of threads to use (default: 1)"
    echo "  --samples    Path to samples CSV file"
    exit 1
}

download_bioproject_reads() {
    local threads=$1
    local samples=$2

    while IFS=',' read -r BioSample BioProject sample library_ID short_R1_filename short_R2_filename long_filename; do

        if [ ! -f "raw_sequencing_reads/${long_filename}" ]; then
            echo "Downloading long reads for $sample"
            fasterq-dump ${long_filename%.fastq.gz} --outdir raw_sequencing_reads --threads $threads
            gzip "raw_sequencing_reads/${long_filename%.gz}"
        else
            echo "Skipping long reads for $sample, already downloaded"
        fi

        if [ ! -f "raw_sequencing_reads/${short_R1_filename}" ] || \
           [ ! -f "raw_sequencing_reads/${short_R2_filename}" ]; then
            echo "Downloading short reads for $sample"
            fasterq-dump --split-files ${short_R1_filename%_1.fastq.gz} --outdir raw_sequencing_reads --threads $threads
            gzip "raw_sequencing_reads/${short_R1_filename%.gz}"
            gzip "raw_sequencing_reads/${short_R2_filename%.gz}"
        else
            echo "Skipping short reads for $sample, already downloaded"
        fi

    done < <(tail -n +2 "$samples")
}

# Defaults
threads=1
samples=""

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --threads) threads="$2"; shift ;;
        --samples) samples="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
    shift
done

# Validate
if [ -z "$samples" ]; then
    echo "Error: --samples is required"
    usage
fi

download_bioproject_reads "$threads" "$samples"