#!/bin/sh

# set your environment here, e.g.:
module load git/2.30.1 python/3.7 perl/5.24.3 snakemake/5.24.1 singularity/3.7.4
module load samtools/1.11

# edit submission command as desired:
cmd="sbatch -t 200:00:00 --export=ALL --mem=32g -p norm -o ${PWD}/SV_wrapper.sh.o%j --wrap='${PWD}/SV_wrapper.sh ${PWD}/config.yaml' "
echo "Command run: $cmd"
eval $cmd