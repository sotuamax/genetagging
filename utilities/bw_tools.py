import pyBigWig 
import pandas as pd 

def bw_value(bw:str, bed:pd.DataFrame, step:int):
    """
    Input: 
    bw: bigwig input file 
    bed: bed df to indicate regions to scan for score
    step: step size used to scan a region 
    Return: 
    mean score over bed regions with a stepsize as fixed value.  
    """
    bw_handle = pyBigWig.open(bw, "r")
    # the default stats is the average value over a range 
    bw_score = [[bw_handle.stats(row.chrom, i, i+step) if i+step < row.end else bw_handle.stats(row.chrom, i, row.end) for i in range(row.start, row.end, step)] for row in bed.itertuples()]
    return bw_score

