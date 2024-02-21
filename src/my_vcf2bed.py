'''
Benchmark the MoCCA-SV results.
Install latest truvari from git clone 
python setup.py install
Finished processing dependencies for Truvari==2.0.0.dev0

So the INSERTION start position is started after the VCF position, and  inserted between two bed posision.
And the Deletion is also starts after the refered position and it dependes on whether ALT is empty or not.

Note SVLEN is the length of ALT minus the length of REF, so negative value suggests DELETION
REF: TCCCCCATGGCTCTTGGCCGTTGGGGCCCAGTTGGCCGCAGAGCCTTTGGCCTTGAGCACTTCCCCAG: len=68
ALT: T
SVLEN=-67
VCF: start 
bed: 22      19947321        19947389        HG3_Ill_manta_4152      DEL     67

REF: G
ALT: GGTCTGCTCACTTGTGGACTGGACACAATTCCTTCTAGGCTGCCGGGGGGAGGATGACAA: len=60
SVLEN=59
VCF: start
Bed: 22      25832627        25832628        HG3_Ill_GATKHC_30298    INS     59
'''
from truvari.utils import *
from MoCCASV.benchmark import *
from MoCCASV.IncludedRegs import IncludedRegs
import argparse
import logging
import os
import pysam



def parse_args():
    """
    Pull the command line parameters
    Refer to: https://github.com/spiralgenetics/truvari/blob/develop/truvari/bench.py
    """
    
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

     
    parser.add_argument("-i", "--infile", type=str, required=True,
                        help="VCF")
    
    args = parser.parse_args()
    # if args.pctsim != 0 and not args.reference:
    #     parser.error("--reference is required when --pctsim is set")

    return args



if __name__ == '__main__':
    
    args=parse_args()
    
    setup_logging(False, LogFileStderr(os.path.join("./", "log.txt")))

    

    base_vcf_fn = args.infile
    logging.info("Process the base vcf file now: %s.", base_vcf_fn)
    base_vcf = pysam.VariantFile(base_vcf_fn) 
    
    base_cnt = 0
    for entry in base_vcf.fetch():
        
        base_type = sv_type(entry)
        base_size = entry_size(entry)
        base_cnt+=1

        if entry.id=='.' or entry.id is None:
            entry.id = "UID_"+str(base_cnt)

        # 1       25204183        .       T       <DUP>   100     PASS    SOMATIC;SVTYPE=DUP;END=25210854;SVLEN=6671      GT      ./.
        # python my_vcf2bed.py -i ../../bench_ws/synthetic/is1.vcf  | grep -P "DUP|INV|DEL" > my_is1.bed

        base_stop=entry.stop
        if entry.stop == entry.start+1 and sv_type != 'INS' and base_size >0 :
            # there is something wrong in parse IS1 data 
            base_stop = entry.start + base_size

        print('\t'.join(map(str, [entry.chrom, entry.start, base_stop, entry.id, base_type,base_size])))

        

    # print(base_cnt)
