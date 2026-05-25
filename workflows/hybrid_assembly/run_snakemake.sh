#!/bin/bash
#$ -S /bin/bash
#$ -pe def_slot 3
#$ -l s_vmem=3G
#$ -l h_rt=72:00:00
#$ -N snakemake_master
#$ -cwd
#$ -o logs/snakemake_master.out

source /home/cgcole/miniforge3/etc/profile.d/conda.sh
conda activate bacterial_genome_pipeline

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_MAIN_FREE=1

snakemake \
    --snakefile workflows/hybrid_assembly/Snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -S /bin/bash -cwd -v PATH=/home/cgcole/miniforge3/bin:/home/cgcole/miniforge3/condabin:\$PATH,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1 -pe def_slot {threads} -l s_vmem={resources.mem_gb}G -l h_rt={resources.runtime} -o logs/{rule}.{jobid}.out -j y" \
    --jobs 10 \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second 1 \
    --latency-wait 60 \
    --use-conda \
    --conda-base-path /home/cgcole/miniforge3 \
    --keep-going \
    --rerun-incomplete \
    > logs/snakemake_hybrid_assembly.log 2>&1 