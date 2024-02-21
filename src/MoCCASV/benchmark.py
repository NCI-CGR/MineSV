
import pandas as pd
from intervaltree import IntervalTree, Interval
from collections import defaultdict, namedtuple
import logging
import pysam
import sys
from MoCCASV.IncludedRegs import IncludedRegs
import re

MATCHRES = namedtuple("matchresult", ("base", "comp", "ovl_pct"))
SVMAT=re.compile("^<(?P<SV>\w+)>$")

def sv_type(entry):
    if "SVTYPE" in entry.info:
        return entry.info["SVTYPE"]
    
    m = SVMAT.match(entry.alts[0])
    if m:
        return m.group('SV')
    
    return None
    

def filter_base_entry(entry, included_regs, sizemin=None, sizemax=None, passonly=False, SV=['INV', 'DUP', 'DEL', 'INS'], SV_ALT= None):
    ''' 
    I may define excluded_regs instead so as to drop any entry overlaps with the excluded regions
    '''

    if passonly and "PASS" not in entry.filter:
        return False

    if SV_ALT is None:
        SV_ALT = map(lambda x: '<'+x+'>', SV)

    is_sv=False

    if entry.alts[0] in SV_ALT or  ("SVTYPE" in entry.info and entry.info["SVTYPE"] in SV):
        is_sv = True

    if not is_sv:
        return False

    if sv_type(entry) in ["INS", "DEL"]:
        # check size only for INDEL
        size = entry_size(entry)

        if sizemin is not None and size<sizemin:
            return False

        if sizemin is not None and size>sizemax:
            return False

    # the base entry should be completely enveloped by the included regions
    # it is a set
    return included_regs.is_included(entry)





def filter_comp_bed(bed_row, included_regs, sizemin=None, sizemax=None):
    chrom, start, stop, id = bed_row

    size = stop-start

    # skip size here as there is no SV type
    # if sizemin is not None and size<sizemin:
    #     return False

    if sizemin is not None and size>sizemax:
        return False

    return included_regs.is_included(bed_row) 

def load_comp_intra(bed_fn, included_regs, sizemin, sizemax, one_based):
    logging.info("Process the intra bed file: %s", bed_fn)
    bed_input = pd.read_csv(bed_fn, header=None, sep="\t", names=['chrom', "start", "end", "id"],dtype={"chrom": "str", "start":"int64", "end":"int64", "id":"int64"})
    
    if one_based:
        bed_input['start'] = bed_input['start'] -1

    
    # print(bed_input.dtypes)
    is_selected = bed_input.apply(lambda x: filter_comp_bed(x, included_regs, sizemin, sizemax) , axis=1)

    comp_cnt = sum(is_selected)
    logging.info("Find %d qualified intra SV predictions.", comp_cnt)
    bed_fil = bed_input[is_selected]

    # make the interval tree for comp bed file
    comp_intra = defaultdict(IntervalTree)
    bed_fil.apply( lambda x: comp_intra[x[0]].addi(x[1],x[2],x[3]), axis=1)

    return comp_cnt, comp_intra

def fetch_coords(lookup, entry, dist=500):
    ''' 
     
    Get the minimum/maximum fetch coordinates to find all variants within dist of variant
    
    similar to fetch_coords defined in truvari
    see https://github.com/spiralgenetics/truvari/blob/e84804b20214775ba625819568e74a1920c99a06/truvari/comparisons.py#L159
    '''
    start, end = entry.start, entry.stop
    start -= dist
    end += dist
    # Membership queries are fastest O(1)
    if not lookup[entry.chrom].overlaps(start, end):
        return None, None

    cand_intervals = lookup[entry.chrom].overlap(start, end)
    s_ret = min([x.begin for x in cand_intervals if x.overlaps(start, end)])
    e_ret = max([x.end for x in cand_intervals if x.overlaps(start, end)])
    return s_ret, e_ret

def reciprocal_overlap(astart, aend, bstart, bend):
    """
    creates a reciprocal overlap rule for matching two entries. Returns a method that can be used as a match operator
    """
    ovl_start = max(astart, bstart)
    ovl_end = min(aend, bend)
    if ovl_start < ovl_end:  # Otherwise, they're not overlapping
        ovl_pct = float(ovl_end - ovl_start) / max(aend - astart, bend - bstart)
    else:
        ovl_pct = 0
    return ovl_pct

def get_match(entry, interval, pad=100):
    '''
    compare the base entry and the comp call
    return None if there is no match 
    and the match information, metrics otherwise

    similar to match_calls() in truvari
    '''
    if not interval.overlaps(entry.start-pad, entry.stop+pad):
        return None

    #print(dir(interval))
    # ovs = interval.overlap_size(entry.start, entry.stop)
    ovl_pct = reciprocal_overlap(entry.start, entry.stop, interval.begin, interval.end)
    return MATCHRES(entry, interval, ovl_pct)


def entry_size(entry):
    """
    https://github.com/spiralgenetics/truvari/blob/ced8a8a1ba8eb16c334de625ad1c902c5f264b8c/truvari/comparisons.py
    returns the size of the variant.
    How size is determined:
    - Starts by trying to use INFO/SVLEN
    - If SVLEN is unavailable and ALT field is an SV (e.g. <INS>, <DEL>, etc),
      use abs(vcf.start - vcf.end). The INFO/END tag needs to be available,
      especially for INS.
    - Otherwise, return the size difference of the sequence resolved call using
      abs(len(vcf.REF) - len(str(vcf.ALT[0])))
    """
    if "SVLEN" in entry.info:
        if type(entry.info["SVLEN"]) in [list, tuple]:
            size = abs(entry.info["SVLEN"][0])
        else:
            size = abs(entry.info["SVLEN"])
    elif entry.alts[0].count("<"):
        start, end = entry_boundaries(entry)
        size = end - start
    else:
        size = abs(len(entry.ref) - len(entry.alts[0]))
    return size