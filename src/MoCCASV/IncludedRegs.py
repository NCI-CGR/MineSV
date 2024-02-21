'''
Modified from https://github.com/spiralgenetics/truvari/blob/develop/truvari/genome_tree.py
'''
import logging 
import pandas as pd
from intervaltree import IntervalTree, Interval
from collections import defaultdict, namedtuple

class IncludedRegs:

    def __init__(self, ref_filename, includebed, chrs= list(map(str, range(1,23)))+['X','Y','MT'] ):
        self.regs=get_included_regs(ref_filename, includebed, chrs)

    
    def is_included(self, entry):
        '''
        Tell whether the bed entry is included in the regs or not
        '''

        entryClass = entry.__class__.__name__
        chrom, start, end = None, None, None

        if entryClass == "VariantRecord":
            chrom, astart, aend = entry.chrom, entry.start, entry.stop
        else: 
            chrom, astart, aend = entry[0:3]

        if chrom not in self.regs: 
            return False

        
        overlaps = self.regs[chrom].overlaps(astart) and self.regs[chrom].overlaps(aend)

        if astart == aend:
            return overlaps

        return overlaps and len(self.regs[chrom].overlap(astart, aend)) == 1

def get_included_regs(ref_filename, includebed, chrs= list(map(str, range(1,23)))+['X','Y','MT']):
    fai_fn = ref_filename+'.fai'

    logging.info("Open the referene sequence index: %s", fai_fn)
    fai_file=open(fai_fn)

    fai_dat = pd.read_csv(fai_fn, sep='\t',header=None)
    #fai_dat[fai_dat[0].isin(chrs)]

    included_regs = defaultdict(IntervalTree)
    fai_dat[fai_dat[0].isin(chrs)].apply( lambda x: included_regs[x[0]].addi(0,x[1]), axis=1)

    if includebed is not None:
        # include region should be the intersection of fai and the include bed file
        bed_regs = bed2tree(includebed, chrs)
        for x in included_regs.keys():
            included_regs[x] = overlap_intervals(included_regs[x], bed_regs[x])

    return included_regs

# This is from https://github.com/chaimleib/intervaltree/commit/3ba0ff3fa2f83643ca71ee5ca70be4d70c59a5cc 
def overlap_intervals(tree1,tree2):
    """
    Returns a new IntervalTree consisting of intervals representing the
    regions overlapped by at least one interval in both of self and other.
    """
    if tree2 is None or tree1 is None:
        return None

    splits = (tree1 | tree2)
    splits.split_overlaps()
    self_int_other =  IntervalTree(filter(lambda r: tree1.overlaps(r) and tree2.overlaps(r), splits))
    self_int_other.merge_overlaps()
    return self_int_other

def bed2tree(bed_fn, chrs):
    logging.info("Create the interval tree for  %s.", bed_fn)
    _dat = pd.read_csv(bed_fn, sep='\t',header=None, dtype={0: "str", 1:"int64", 2:"int64"})
    
    _regs = defaultdict(IntervalTree)
    _dat[_dat[0].isin(chrs)].apply( lambda x: _regs[x[0]].addi(x[1],x[2]), axis=1)
    return _regs