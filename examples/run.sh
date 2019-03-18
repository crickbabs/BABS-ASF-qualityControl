#!/bin/sh

module load Anaconda2

source activate cotidianus

command="sbatch"
command="$command --cpus-per-task 1"
command="$command --nodes 1"
command="$command --ntasks 1" 
command="$command --output slurm.out"
command="$command --error slurm.err"
command="$command --mem-per-cpu 6G"

snakemake --jobs 999 --cluster "$command" --rerun-incomplete

source deactivate

