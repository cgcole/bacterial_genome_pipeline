# Bacterial Genome Assembly Pipeline

This pipeline can be used to hybrid assemble short and long read whole genome sequencing using [Autocycler](https://github.com/rrwick/Autocycler/wiki). Currently, you can assemble internal raw sequencing files, or you can specify BioProjects to download and assemble publicly available raw whole genome sequencing reads.

## Requirements

- conda or mamba

## Installation

1. Clone the repository:
```bash
    git clone https://github.com/cgcole/bacterial_genome_pipeline.git
    cd bacterial_genome_pipeline
```

2. Create and activate the conda environment:
```bash
    conda env create -f envs/environment.yaml
    conda activate bacterial_genome_pipeline
```

3. Download the Plassembler database for Autocycler:
```bash
    bash scripts/plassembler_db_download.sh
```

## Usage

### Downloading metadata for BioProjects of interest
```bash
python scripts/bioproject_metadata_parser.py PRJNAXXXXXX
```

### Downloading raw sequencing reads for BioProjects
```bash
bash scripts/download_bioproject_reads.sh --samples samples/hybrid_assembly/bioproject_hybrid_samples.csv --threads 1
bash scripts/download_bioproject_reads.sh --samples samples/hybrid_assembly/bioproject_hybrid_samples.csv --threads 4
```

### Downloading raw sequencing reads for BioProjects on a cluster
```bash
mkdir -p logs raw_sequencing_reads
NSAMPLES=$(( $(wc -l < samples.csv) - 1 ))
qsub -t 1-${NSAMPLES} download_reads.sh
```

### Assembling genomes with Snakemake
```bash
snakemake --cores x
```

### Assembling genomes on a cluster
```bash
snakemake \
    --executor slurm \
    --default-resources slurm_partition=general mem_mb=16000 runtime=120 \
    --jobs 100
```

## Running the pipeline

### Parameters you may want to adjust:

| Parameter | What it does | Default |
|-----------|-------------|---------|
| `--jobs` | Max number of jobs running simultaneously | 50 |
| `mem_gb` | Memory per slot (GB) | 4 |
| `runtime` | Max time per job (HH:MM:SS) | 12:00:00 |
| `{threads}` | Threads per job, set in config.yaml | 4 |

### Real run on the HPC cluster:
snakemake \
  --snakefile workflows/hybrid_assembly/Snakefile \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "qsub -S /bin/bash -pe def_slot {threads} -l s_vmem={resources.mem_gb}G -l h_rt={resources.runtime}:00 -o logs/{rule}.out -j y" \
  --jobs 50 \
  --latency-wait 60 \
  -n  # remove this -n for real run

### Real run locally (laptop/desktop):
snakemake \
  --snakefile workflow/Snakefile \
  --cores 8 \
  --use-conda