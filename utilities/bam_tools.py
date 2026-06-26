import pysam 
import subprocess
import os 
import pandas as pd 
import bioframe as bf 



def bed2bam(bed, bam, new_bam, n = 4):
    """
    bam servers as a header template; 
    bed provide read coordinate information; 
    write into a new bam file which stores read information w/ cooridnates from bed and other filled by pseudo-information. 
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    bed_df = bf.read_table(bed, schema = "bed6")
    with pysam.AlignmentFile(new_bam, "wb", template = bam_handle, threads = n) as bam_write:
        for row in bed_df.itertuples():
            # only use read as a placeholder, other information is not meanful 
            read = pysam.AlignedSegment()
            read.query_name = row.name
            read.query_sequence = "N"*(row.end-row.start)  # must match CIGAR length
            read.flag = 0 if row.strand == "+" else 16 # mapped, forward strand
            read.reference_id = bam_write.get_tid(row.chrom)
            read.reference_start = row.start  # 0-based start
            read.mapping_quality = row.score
            read.cigar = [(0, row.end - row.start)]  # 0 = 'M' (match)
            read.next_reference_id = -1
            read.next_reference_start = -1
            read.template_length = 0
            read.query_qualities = pysam.qualitystring_to_array("F"*(row.end - row.start))
            bam_write.write(read)

def flagstats(bam, n):
    name = bam.replace(".bam", "")
    subprocess.call(f"samtools flagstats {bam} -@ {n} -O tsv > {name}.flagstat", shell = True)

def examine_sort(bam):
    """
    Examine BAM as name or coordinates sorted
    """
    bamfile = pysam.AlignmentFile(bam, "rb")
    # Access the header dictionary
    header = bamfile.header
    # Check sorting order from @HD line
    sort_order = header.get('HD').get('SO')
    # standard output: queryname (sort by name); coordinate (sort by coordinate)
    return sort_order 

def nsort_clean(bam:str, new_bam:str, n = 1):
    """
    When QC applied on each read, read and its mate may not be retained at the same time. This func is to clean BAM file for reads that forms in pair. 
    Input: 
    bam: name sorted bam file 
    new_bam: new bam file name to write into
    n: threads to use
    Return: 
    write into new bam file with read and its mate sorted by name. 
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    for read in bam_handle:
        read0 = read 
        break
    with pysam.AlignmentFile(new_bam, "wb", template = bam_handle, threads = n) as newbam:
        for read in bam_handle:
            if read.query_name == read0.query_name:
                newbam.write(read0)
                newbam.write(read)
            read0 = read 
    
def cigar_extract(read_df):
    """
    Extract from cigarstring for cigar label and its matched length
    Input: 
    read_df (dataframe): read bed extracted from BAM file 
    Return: 
    read with cigar informaiton in a dataframe (index:read_name, cigar_char, cigar_len)
    """
    soft_chars = read_df["cigar"].str.findall("[M|N|S]+")
    cigar_digit = read_df["cigar"].str.findall("[0-9]+")
    cigar_df = pd.concat([soft_chars, cigar_digit], axis = 1)
    cigar_df.columns = ["cigar_char", "cigar_len"]
    #cigar_df["cigar_len"] = cigar_df["cigar_len"].astype(int)
    cigar_df.index = read_df.index
    return cigar_df

def bam2bed(bam, label = "cDNA", read_seq = False, chrom = None):
    """
    Transform bam into bed format
    Input: 
    bam: bam file for RNA-seq data
    label: bam file type (cDNA/DNA)
    Output: 
    read bed file (chrom, start, end, name, score, strand, cigar, read-seq as optional)
    """
    bam_handle = pysam.AlignmentFile(bam, "rb")
    strandness = {True:"+", False:"-"}
    # chrom, start, end, name, score, strandness, cigar
    # read_df = pd.DataFrame(read_df, columns = ["read", "seq", 'strand', "cigar"])
    if label.upper() in ["CDNA", "RNA"]:
        # for RNA-seq alignment, ("NH", 1) to get unique alignment (no quality filter, MAPQ see score)
        if read_seq:
            read_bed = pd.DataFrame([[read.reference_name, read.reference_start,read.reference_end, read.query_name, strandness[read.is_forward], read.cigarstring, read.query_sequence] for read in bam_handle.fetch(chrom) if ("NH", 1) in read.tags])
        else:
            read_bed = pd.DataFrame([[read.reference_name, read.reference_start,read.reference_end, read.query_name, strandness[read.is_forward]] for read in bam_handle.fetch(chrom) if read.get_tag("NH") == 1])
    if label.upper() == "DNA":
        # for DNA-seq alignment, filtering on primary alignment (no quality filter, MAPQ see score)
        if read_seq:
            read_bed = pd.DataFrame([[read.reference_name,read.reference_start,read.reference_end,read.query_name, strandness[read.is_forward], read.cigarstring, read.query_sequence] for read in bam_handle.fetch() if (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped)])
        else:
            read_bed = pd.DataFrame([[read.reference_name,read.reference_start,read.reference_end,read.query_name, strandness[read.is_forward], read.cigarstring] for read in bam_handle.fetch() if (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped)])
    if read_seq and len(read_bed) > 0:
        read_bed.columns = ["chrom", "start", "end", "name", "strand", "cigar", "seq"]
        return read_bed
    if not read_seq and len(read_bed) > 0:
        read_bed.columns = ["chrom", "start", "end", "name", "strand"]
        return read_bed
    if len(read_bed) == 0:
        return None

def sc_markdup(bam, output, n):
    """
    Perform markdup for single cell BAM data 
    output: 
    write new BAM (de-duplicated) to output; 
    write fragment file (for each of the fragment, label its frquency); 
    duplication rate; 
    """
    strandness = {True:"+", False:"-"}
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    i = 0; unique_i = 0
    with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = n) as newbam:
        for chr in bam_header(bam, seq = True):
            print(f"Processing {chr} ...")
            bam_tmp = [(read.query_name.split(":")[-1], read.reference_start, read.reference_end, read.query_name, strandness[read.is_forward], read) for read in bam_handle.fetch(chr)]
            if len(bam_tmp) > 0:
                bam_tmp = pd.DataFrame(bam_tmp, columns = ["index", "start", "end", "name", "strand", "read"])
                bam_tmp["label"] = bam_tmp.groupby(["index", "start", "end", "strand"])["name"].rank(method = "first").astype(int)
                bam_tmp["chrom"] = chr; bam_tmp["n"] = bam_tmp.groupby(["index", "start", "end", "strand"]).transform("size")
                #bam_tmp[["chrom", "start", "end", "strand", "index", "n"]].drop_duplicates(keep = "first").to_csv(output.replace(".bam", ".fragment"), sep = "\t", header = False, index = False, mode = "a")
                bam_pick = bam_tmp.query("label == 1")
                for row in bam_pick.itertuples():
                    newbam.write(row.read)
                i += len(bam_tmp)
                unique_i += len(bam_pick)
    subprocess.call(f"samtools index -@ {n} {output}", shell = True)
    df_dup = pd.DataFrame([os.path.basename(bam), i, unique_i, 1-unique_i/i], index = ["LIBRARY", "UNPAIRED_READS_EXAMINED", "UNIQUE_READS", "PERCENT_DUPLICATION"], columns = ["stat"])
    df_dup.to_csv(output.replace(".bam", ".stat"), sep = "\t", header = True, index = True)

def bam2bedpe(bam, name):
    """
    Transform bam into bedpe format using bedtools 
    Input: 
    bam: bam file 
    Output: 
    bedpe file 
    """
    # chrom, start, end, name, score, strandness, cigar
    # read_df = pd.DataFrame(read_df, columns = ["read", "seq", 'strand', "cigar"])
    bedpe_command = f"bedtools bamtobed -i {bam} -bedpe > {name}"
    subprocess.call(bedpe_command, shell = True)

def bam2bw(bam, bw, n, smooth = 20, norm = False):
    if norm:
        bw_command = f"bamCoverage -b {bam} -o {bw} -of bigwig -p {n} --normalizeUsing CPM --smoothLength {smooth}"
    else:
        bw_command = f"bamCoverage -b {bam} -o {bw} -of bigwig -p {n} --smoothLength {smooth}"
    subprocess.call(bw_command, shell=True)

def bam_count(bam, chr = None, start = None, end = None, n = 1):
    """
    Count number of reads within specified region of given bam file
    Input: 
    bam: bam file 
    chr: chromosome region 
    start: chromosome start region 
    end: chromosome end region 
    n: thread to process bam file 
    Output: 
    number of reads within the specified region of bam file
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    read_set = {read.query_name for read in bam_handle.fetch(chr, start, end)}
    # for read in bam_handle.fetch()
    return len(read_set)

def bam_count_fast(bam, chr = None, start = None, end = None, n = 1):
    """
    For PE reads, count total number of reads and then divide it by 2 to get estimated read number
    # Note: it is not accurate, as some paired-reads only have 1-end of it pass quality control. 
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    read_count = bam_handle.count(chr, start,end)
    return read_count//2

def bam_merge(bams:list, name, thread):
    file_combined = " ".join(bams)
    merge_command = f"samtools merge -r -@ {thread} -o - {file_combined} | samtools sort -@ {thread} - -o {name}.bam && samtools index -@ {thread} {name}.bam"
    if not os.path.exists(name + ".bam"):
        print(merge_command)
        subprocess.call(merge_command, shell = True)

def bam_header(bam, seq = False):
    """
    Give BAM file, to return BAM header / BAM involved sequences 
    """
    bam_handle = pysam.AlignmentFile(bam, "rb")
    bam_header = str(bam_handle.header)
    if seq:
        bam_seq = [s.split("\t")[1].split(":")[-1] for s in str(bam_handle.header).split("\n") if s.startswith("@SQ")]
        return bam_seq 
    else:
        return bam_header 

def bed2saf(bed):
    """
    Transfer bed file input saf format (see featureCount)
    Input: 
    bed: bed file 
    Output:
    None
    """
    import bioframe as bf 
    name = os.path.basename(bed).replace("bed", "saf")
    bed_df = bf.read_table(bed, schema = "bed5")
    saf = bed_df[["name", "chrom", "start", "end"]]
    saf["strand"] = "."
    saf.to_csv(name, sep = "\t", header = False, index = False)

def bamfilter(bam, output, select_read = None, threads = 1, clip_check = False, two_end_clip_check = False, chrom = False):
    """
    Rewrite new bam file by given selected reads name
    """
    # strandness = {True:"+", False:"-"}
    if not output.endswith(".bam"):
        output += ".bam"
    if chrom:
        import bioframe as bf 
        refseq = bf.assembly_info("hg38")
        ref_list = refseq.seqinfo["name"].tolist()
        ref_list.remove("chrY"); ref_list.remove("chrM")
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = threads)
    if select_read is not None:
        with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
            for read in bam_handle.fetch():
                if read.query_name in select_read:
                    newbam.write(read)
    else:
        if clip_check and chrom:
            import re
            # need to examine the 5' end read for any soft-clip 
            with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
                for read in bam_handle.fetch():
                    if read.reference_name in ref_list and read.mapping_quality > 30 and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
                        if "S" in read.cigarstring:
                            if (re.findall(r'[A-Z=]', read.cigarstring)[0] == "S" and read.is_forward) or (re.findall(r'[A-Z=]', read.cigarstring)[-1] == "S" and read.is_reverse):
                                pass 
                            else:
                                newbam.write(read)
                        else:
                            newbam.write(read)
        elif clip_check and not chrom:
            print("Filter on 5' end (no clip) ...")
            import re
            # need to examine the 5' end read for any soft-clip 
            with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
                for read in bam_handle.fetch():
                    if read.mapping_quality > 30 and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
                        if "S" in read.cigarstring:
                            if (re.findall(r'[A-Z=]', read.cigarstring)[0] == "S" and read.is_forward) or (re.findall(r'[A-Z=]', read.cigarstring)[-1] == "S" and read.is_reverse):
                                pass 
                            else:
                                newbam.write(read)
                        else:
                            newbam.write(read)
        elif two_end_clip_check and chrom:
            with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
                for read in bam_handle.fetch():
                    if read.reference_name in ref_list and read.mapping_quality > 30 and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
                        if (not "S" in read.cigarstring):
                            newbam.write(read)
        elif two_end_clip_check and not chrom:
            with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
                for read in bam_handle.fetch():
                    if read.mapping_quality > 30 and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
                        if ("S" not in read.cigarstring):
                            newbam.write(read)
        elif not clip_check and chrom: 
            with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
                for read in bam_handle.fetch():
                    if read.reference_name in ref_list and read.mapping_quality > 30 and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
                        newbam.write(read)
        else:
            with pysam.AlignmentFile(output, "wb", template = bam_handle, threads = threads) as newbam:
                for read in bam_handle.fetch():
                    if read.mapping_quality > 30 and (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped):
                        newbam.write(read)
    try:
        subprocess.call(f"samtools index -@ {threads} {output}", shell = True)
    except Exception as e:
        print(e)

def feature_count(bam, feature_df):
    import bioframe as bf 
    bam_handle = pysam.AlignmentFile(bam, "rb")
    read_list = [(read.reference_id, read.reference_start, read.reference_end, read.query_name) for read in bam_handle.fetch()]
    read_df = pd.DataFrame(read_list, columns = ["chrom", "start", "end", "name"])
    read_overlap_feature = bf.overlap(read_df, feature_df) 
    # feature_df: exon, intro, intergenic region 
    return read_overlap_feature 

