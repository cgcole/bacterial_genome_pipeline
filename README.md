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
bash scripts/download_bioproject_reads.sh --samples samples/hybrid_assembly/bioproject_hybrid_samples.csv --threads 4

```

### Downloading raw sequencing reads for BioProjects on a cluster
```bash
mkdir -p logs raw_sequencing_reads
NSAMPLES=$(( $(wc -l < samples.csv) - 1 ))
qsub -t 1-${NSAMPLES} scripts/download_bioproject_reads.sh

mkdir -p logs raw_sequencing_reads
NSAMPLES=$(( $(wc -l < samples/short_assembly/bioproject_short_samples.csv) - 1 ))
qsub -t 1-${NSAMPLES} scripts/download_bioproject_reads.sh
```

## Running the pipeline

### Real run on the HPC cluster:
snakemake \
  --snakefile workflows/hybrid_assembly/Snakemake \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "qsub -S /bin/bash -cwd -v PATH=/home/cgcole/miniforge3/bin:/home/cgcole/miniforge3/condabin:$PATH,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1 -pe def_slot {threads} -l s_vmem={resources.mem_gb}G -l h_rt={resources.runtime} -o logs/{rule}.{jobid}.out -j y" \
  --jobs 10 \
  --latency-wait 60 \
  --use-conda \
  --conda-base-path /home/cgcole/miniforge3 \
  -n  # remove this -n for real run

  ### Parameters you may want to adjust:

| Parameter | What it does | Default |
|-----------|-------------|---------|
| `--jobs` | Max number of jobs running simultaneously | 50 |
| `mem_gb` | Memory per slot (GB) | 4 |
| `runtime` | Max time per job (HH:MM:SS) | 12:00:00 |
| `{threads}` | Threads per job, set in config.yaml | 4 |

### Real run locally (laptop/desktop):
snakemake \
  --snakefile workflow/Snakefile \
  --cores 8 \
  --use-conda \
  -n  # remove this -n for real run










nohup snakemake \
    --snakefile workflows/hybrid_assembly/Snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -S /bin/bash -cwd -v PATH=/home/cgcole/miniforge3/bin:/home/cgcole/miniforge3/condabin:\$PATH,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1 -pe def_slot {threads} -l s_vmem={resources.mem_gb}G -l h_rt={resources.runtime} -o logs/{rule}.{jobid}.out -j y" \
    --jobs 10 \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second .1 \
    --latency-wait 60 \
    --use-conda \
    --conda-base-path /home/cgcole/miniforge3 \
    --keep-going \
    --rerun-incomplete \
    > logs/snakemake_hybrid_assembly.log 2>&1 &



















nohup snakemake \
    --snakefile workflows/short_assembly/Snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -S /bin/bash -cwd -v PATH=/home/cgcole/miniforge3/bin:/home/cgcole/miniforge3/condabin:\$PATH,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1 -pe def_slot {threads} -l s_vmem={resources.mem_gb}G -l h_rt={resources.runtime} -o logs/{rule}.{jobid}.out -j y" \
    --jobs 10 \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second .1 \
    --latency-wait 60 \
    --use-conda \
    --conda-base-path /home/cgcole/miniforge3 \
    --keep-going \
    > logs/snakemake_short_assembly.log 2>&1 &









nohup snakemake \
    --snakefile workflows/download_bioproject_reads/Snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -S /bin/bash -cwd -v PATH=/home/cgcole/miniforge3/bin:/home/cgcole/miniforge3/condabin:\$PATH,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1 -pe def_slot {threads} -l s_vmem={resources.mem_gb}G -l h_rt={resources.runtime} -o logs/{rule}.{jobid}.out -j y" \
    --jobs 10 \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second .1 \
    --latency-wait 60 \
    --use-conda \
    --conda-base-path /home/cgcole/miniforge3 \
    --keep-going \
    --rerun-incomplete \
    > logs/snakemake_download_bioproject_reads.log 2>&1 &