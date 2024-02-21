'''
Benchmark the MoCCA-SV results.

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
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

     
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline truth-set VCF calls")
    parser.add_argument("-c", "--comp", type=str, required=True,
                        help="Comparison set of calls in bed file")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Indexed fasta used to call variants")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")
    parser.add_argument("--onebased", action="store_true", default=False,
                        help="Intra bed is 1-based")    
    # parser.add_argument("--prog", action="store_true",
    #                    help="Turn on progress monitoring")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=restricted_float, default=0.7,
                        help="Minimum pct reciprocal overlap (%(default)s) for SV events")


    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizefilt", type=int, default=30,
                        help="Minimum variant size to load into IntervalTree (%(default)s)")
    filteg.add_argument("--sizemax", type=int, default=50000,
                        help="Maximum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("--includebed", type=str, default=None,
                        help="Bed file of regions in the genome to include only calls overlapping")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")
    filteg.add_argument("--multimatch", action="store_true", default=False,
                        help=("Allow base calls to match multiple comparison calls, and vice versa. "
                              "Output vcfs will have redundant entries. (%(default)s)"))

    args = parser.parse_args()
    # if args.pctsim != 0 and not args.reference:
    #     parser.error("--reference is required when --pctsim is set")

    return args



if __name__ == '__main__':
    
    args=parse_args()
    
    setup_logging(args.debug, LogFileStderr(os.path.join("./", "log.txt")))

    # -f ../../refGenomes/Homo_sapiens_assembly19.fasta
    # -b "../data/dream_modified.vcf.gz"
    # -c "../data/IS1_intra.bed"
    # included_reg = get_included_reg(args.reference, args.includebed)
    included_reg = IncludedRegs(args.reference, args.includebed)
    
    intra_cnt, intra_tree=load_comp_intra(args.comp, included_reg, args.sizefilt, args.sizemax, args.onebased)

    base_vcf_fn = args.base
    logging.info("Process the base vcf file now: %s.", base_vcf_fn)
    base_vcf = pysam.VariantFile(base_vcf_fn) 
    
    tp=0
    base_cnt=0

    for entry in base_vcf.fetch():
        if not filter_base_entry(entry, included_reg, args.sizemin, args.sizemax, args.passonly) :
            continue

        base_type = sv_type(entry)
        base_cnt+=1

        

        fetch_start, fetch_end = fetch_coords(intra_tree, entry, args.refdist)
        
        

        if fetch_start is None and fetch_end is None:
            continue

        mats = []
        # now get the matched calls from the given region
        for comp_sv in intra_tree[entry.chrom].overlap(fetch_start, fetch_end):
            mat = get_match(entry, comp_sv, args.refdist)
            
            if mat is None:
                continue

            # print(mat.base.chrom, mat.base.start, mat.base.stop, mat.comp.begin, mat.comp.end, mat.ovl_pct)
            mats.append(mat)

        # select the best match
        if len(mats)>0:
            mats.sort(reverse=True, key=lambda x:x.ovl_pct)
            mat = mats[0]
            if mat.ovl_pct < args.pctovl and base_type=="DEL":
                continue

            logging.info("Match: %s\t%s\t%d\t%d\t%d\t%d\t%7.6f", base_type, mat.base.chrom, mat.base.start, mat.base.stop, mat.comp.begin, mat.comp.end, mat.ovl_pct)
            tp+=1

    print(base_cnt, intra_cnt, tp)
