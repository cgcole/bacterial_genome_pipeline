#!/bin/bash
#$ -l s_vmem=4G  # or whatever resources you need
#$ -l h_rt=12:00:00
#$ -o logs/download.$TASK_ID.out
#$ -e logs/download.$TASK_ID.err
#$ -j y   # or merge stdout and stderr into one file

source /home/cgcole/miniforge3/etc/profile.d/conda.sh
conda activate bacterial_genome_pipeline


usage() {
    echo "Usage: $0 --threads <int> --samples <file> --outdir <dir>"
    echo "  --threads    Threads per job (default: 4)"
    echo "  --samples    Path to samples CSV file"
    echo "  --outdir     Output directory for reads (default: raw_sequencing_reads)"
    exit 1
}

# ── If running as a cluster array task, skip straight to download ──
if [ -n "$SGE_TASK_ID" ]; then

    if [ -z "$SAMPLES" ] || [ -z "$READS_DIR" ] || [ -z "$THREADS" ]; then
        echo "Error: SAMPLES, READS_DIR, and THREADS must be set" >&2
        exit 1
    fi

    IFS=',' read -r BioSample BioProject sample library_ID short_R1 short_R2 long_file \
        < <(sed -n "$((SGE_TASK_ID + 1))p" "$SAMPLES")

    if [ ! -f "${READS_DIR}/${long_file}" ]; then
        echo "Downloading long reads for $sample"
        fasterq-dump "${long_file%.fastq.gz}" \
            --outdir "$READS_DIR" \
            --threads "$THREADS" && \
        pigz -p "$THREADS" "${READS_DIR}/${long_file%.gz}"
    else
        echo "Skipping long reads for $sample, already downloaded"
    fi

    if [ ! -f "${READS_DIR}/${short_R1}" ] || [ ! -f "${READS_DIR}/${short_R2}" ]; then
        echo "Downloading short reads for $sample"
        fasterq-dump --split-files "${short_R1%_1.fastq.gz}" \
            --outdir "$READS_DIR" \
            --threads "$THREADS" && \
        pigz -p "$THREADS" \
            "${READS_DIR}/${short_R1%.gz}" \
            "${READS_DIR}/${short_R2%.gz}"
    else
        echo "Skipping short reads for $sample, already downloaded"
    fi

    exit 0
fi

# ── Otherwise, parse arguments and submit or run locally ──
threads=4
samples=""
outdir="raw_sequencing_reads"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --threads) threads="$2"; shift ;;
        --samples) samples="$2"; shift ;;
        --outdir)  outdir="$2";  shift ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
    shift
done

if [ -z "$samples" ]; then
    echo "Error: --samples is required"
    usage
fi

NSAMPLES=$(( $(wc -l < "$samples") - 1 ))
mkdir -p logs "$outdir"

if command -v qsub &>/dev/null; then
    echo "Submitting $NSAMPLES jobs to cluster..."
    qsub -t 1-${NSAMPLES} \
         -pe def_slot "$threads" \
         -v SAMPLES="$(realpath $samples)",READS_DIR="$(realpath $outdir)",THREADS="$threads" \
         "$0"    # passes this script itself to qsub
else
    echo "No cluster detected, running locally with max 10 parallel jobs..."
    seq 1 $NSAMPLES | parallel -j 10 \
        SAMPLES="$samples" READS_DIR="$outdir" THREADS="$threads" SGE_TASK_ID={} bash "$0"
fi