import pysam 
import pandas as pd 

def read_count(fastq):
    """
    Count read number in fastq 

    Input: 
    fastq: fastq file
    Return:
    read count
    """
    i = 0
    with pysam.FastxFile(fastq) as q_handle: 
        for entry in q_handle:
            i += 1
    return i

def read_len(fastq):
    """
    taste read length in fastq

    Input: 
    fastq: fastq file
    Output:
    read length
    """
    with pysam.FastxFile(fastq) as q_handle:
        first_read = next(q_handle, None)
        if first_read is not None:
            return len(first_read.sequence)

def fq2df(read, quality = False, out = None):
    """
    Transform fastq file into a dataframe, with columns:
    name: read name 
    comment: read description
    seq: read seq 
    quality: (optional) seq quality
    """
    with pysam.FastxFile(read) as q_handle:
        if next(q_handle, None) is None:
            return None # return None if the iterator reaches the end 
    if out is None:
        with pysam.FastxFile(read) as q_handle:
            if quality: 
                read_df = pd.DataFrame([[entry.name, entry.comment, entry.sequence, entry.quality] for entry in q_handle]) 
                read_df.columns = ["name", "comment", "seq", "quality"]
            else:
                read_df = pd.DataFrame([[entry.name, entry.comment, entry.sequence] for entry in q_handle]) 
                read_df.columns = ["name", "comment", "seq"]
            return read_df 
    else:
        with pysam.FastxFile(read) as q_handle, open(out, "w") as tmpw:
            if quality:
                tmpw.write(f"name\tcomment\tseq\tquality\n")
                for entry in q_handle:
                    tmpw.write(f"{entry.name}\t{entry.comment}\t{entry.sequence}\t{entry.quality}\n")
            else:
                tmpw.write(f"name\tcomment\tseq\n")
                for entry in q_handle:
                    tmpw.write(f"{entry.name}\t{entry.comment}\t{entry.sequence}\n")

def read_select(read_df:pd.DataFrame, supporting_pattern:list):
    """
    To select reads match specific pattern at its beginning (5' end).

    Input: 
    read_df: read df with "seq" column
    supporting_pattern: list (in list of strings, all pattern(s) only for 1 length)

    Return: 
    a read df agree with supporting pattern(s)
    """
    read_ddf = read_df.copy()
    read_ddf['p'] = read_ddf["seq"].str.slice(0, len(supporting_pattern[0]))
    if len(supporting_pattern) == 1:
        read_select = read_ddf[read_ddf["p"] == supporting_pattern[0]]
    else:
        read_select = read_ddf[read_ddf['p'].isin(supporting_pattern)]
    return read_select.drop("p", axis = 1)

def write_read(select_df:pd.DataFrame, start:int, output:str, end = None, comment = False):
    """
    write fastq read in a dataframe to fastq file, cut 5' end site when assigned with cutting size (start)
    Input: 
    select_df: df for selected reads 
    start: int (before start site - cut, after start site - keep)
    end: control read size (optional)
    output: output file name 

    Return: 
    None (directly write to file)
    """
    select_ddf = select_df.copy()
    if start != 0:
        #select_ddf["cut"] = select_ddf["seq"].str.slice(0, start)
        # update sequence for write into new file 
        select_ddf["seq"] = select_ddf["seq"].str.slice(start, end)
        select_ddf["quality"] = select_ddf["quality"].str.slice(start, end)
        if output.endswith(".gz"):
            import gzip
            if not comment:
                with gzip.open(output, "wt") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name} {row.comment}\n{row.seq}\n+\n{row.quality}\n")
            else:
                with gzip.open(output, "wt") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name}:{row.comment}\n{row.seq}\n+\n{row.quality}\n")
        else:
            if not comment:
                with open(output, "w") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name} {row.comment}\n{row.seq}\n+\n{row.quality}\n")
            else:
                with open(output, "w") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name}:{row.comment}\n{row.seq}\n+\n{row.quality}\n")
    else:
        if output.endswith(".gz"):
            import gzip
            if not comment:
                with gzip.open(output, "wt") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name} {row.comment}\n{row.seq}\n+\n{row.quality}\n")
            else:
                with gzip.open(output, "wt") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name}:{row.comment}\n{row.seq}\n+\n{row.quality}\n")
        else:
            if not comment:
                with open(output, "w") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name} {row.comment}\n{row.seq}\n+\n{row.quality}\n")
            else:
                with open(output, "w") as fr:
                    for row in select_ddf.itertuples():
                        fr.write(f"@{row.name}:{row.comment}\n{row.seq}\n+\n{row.quality}\n")

def write_head(select_df:pd.DataFrame, start:int):
    """
    Capture head-seq of seq in a dataframe. head-size should be given in "start". 
    """
    select_df = select_df.copy()
    select_df["head"] = select_df["seq"].str.slice(0, start)
    return select_df[["name", "head"]]

def select_regex(read_df:pd.DataFrame, primer_seq:str, regex = True):
    """
    To select reads containing specific primer-seq at its beginnng (5' end).
    Input: 
    read_df: read in a df with "seq" column 
    primer_seq: str

    Return: 
    a df with reads matched to the designated primer
    """
    iupac_codes = {
        'A': r'A',
        'C': r'C',
        'G': r'G',
        'T': r'T',
        'R': r'[A|G]',
        'Y': r'[C|T]',
        'K': r'[G|T]',
        'M': r'[A|C]',
        'B': r'[C|G|T]',
        'D': r'[A|G|T]',
        'H': r'[A|C|T]',
        'V': r'[A|C|G]',
        'W': r'[A|T]',
        'S': r'[C|G]',
        'N': r'[A|C|G|T]',
        'X': r'[A|C|G|T]'
    }
    read_df["primer"] = read_df["seq"].str.slice(0, len(primer_seq))
    primer_seq = primer_seq.upper()
    if regex:
        import re 
        code_pattern = "".join([iupac_codes[bp] for bp in primer_seq])
        read_df_select = read_df[read_df["primer"].str.match(code_pattern, regex = True)].drop(labels = "primer", axis = 1)
    else:
        read_df_select = read_df[read_df["primer"] == primer_seq].drop(labels = "primer", axis = 1)
    return read_df_select
