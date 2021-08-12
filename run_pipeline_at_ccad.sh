#!/bin/sh

# set your environment here, e.g.:
module load  git/2.12.0 python3/3.7.0 sge/2011.11p1 perl/5.32.1 singularity/3.0.1

# edit submission command as desired:
cmd="qsub -q long.q -V -j y -S /bin/sh -o ${PWD} ${PWD}/SV_wrapper.sh ${PWD}/config.yaml"
echo "Command run: $cmd"
eval $cmd
