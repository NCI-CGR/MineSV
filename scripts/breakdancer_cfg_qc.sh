#!/bin/bash

# This script performs the QC checks on the BAMs and config file produced by breakdancer that are described in the tool's documentation.

my_config=$1
#n_bam=$3
out=$2

# Check if the configuration file generated by bam2cfg.pl contains one line for each read group or library in every input BAM file.
# If the first number is smaller than the second one, then this input BAM file may be missing RG or LB information in the header.

echo "Config file QC:" >> $out

# Check if the configuration file contains “NA” as read group.
# If output is greater than 0, then the RG or LB information is missing in at least one of the BAM files.
num3=$(grep -cw readgroup:NA $my_config)
if [ "$num3" -gt "0" ]; then
	echo "FAIL: one of the BAMs is missing RG or LB information." >> $out
else
	echo "PASS: all RG or LB information is present in the bams." >> $out
fi

# Check the coefficient of variation (standard deviation divided by mean) of the insert size for each read group.
# The command prints the coefficient of variation for each read group. Normally, they should be < 0.2 or 0.3.
echo "Coefficient of variation of the insert size for each read group (expected 0.2 or 0.3):" >> $out
perl -ane ' ($mean)=($_=~/mean:(\S+)/);($std)=($_=~/std:(\S+)/);print $std/$mean ."\n" ' $my_config >> $out

# Check the percentage of inter-chromosomal read pairs.
# The above command prints the percentage of the reads with tag 32 for each read group, which corresponds to inter-chromosomal read pairs. This percentage should normally be smaller than 3%.
echo "Percentage of inter-chromosomal read pairs (expected <3%):" >> $out
perl -ane ' ($CTX)=($_=~/32\((\S+?)\)/);print $CTX."\n" ' $my_config >> $out
