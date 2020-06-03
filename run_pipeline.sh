#!/bin/sh

# set your environment here, e.g.:
module load python3/3.7.0 sge perl singularity/3.0.1 git

# edit submission command as desired:
cmd="qsub -q long.q -V -j y -S /bin/sh -o ${PWD} path/to/SV_wrapper.sh ${PWD}/config.yaml"
echo "Command run: $cmd"
eval $cmd
