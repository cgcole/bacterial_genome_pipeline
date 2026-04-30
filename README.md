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


