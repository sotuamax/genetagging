#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program is to clean sequencing data based on the primer given (R1, enriched), 
and assign cell identify based on cell barcode given in R2.
Output file is clean and cell barcode labeled R1 fastq file 
1). start with primer; 
2). valid cell barcode. 
"""

import pysam 
import os 
import sys
import pandas as pd
import collections
import argparse
import logging
import numpy as np
import bioframe as bf
import subprocess
from utilities import * 
from Bio.Seq import reverse_complement
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("-assembly", "--assembly", required = False, help = "reference genome assembly")
    parser.add_argument("-dtype", "--dtype", default = "cDNA", choices = ["cDNA", "DNA"], type = str, help = "data type of the sequencing file used for alignment")
    parser.add_argument("-o", "--output", required = True, help="output directory name")
    parser.add_argument("-gtf", "--gtf", required = False, help = "reference gtf file used for alignment")
    # separate argumentparser for sub-parsers
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    # initiate sub-parser 
    sub_parsers = parser2.add_subparsers(dest = 'command', help = "mode to run")
    # pre-processing
    pre = sub_parsers.add_parser("pre", help = "preprocess fastq sequencing read", parents = [parser], add_help = False)
    pre.add_argument("-fastq", "--fastq", nargs = 2, required = True, help="fastq input file (for paired-end reads, order in R1 R2)")
    pre.add_argument("-read_len", "--read_length", type = int, default = 35, help = "read length used for alignment")
    pre.add_argument("-primer", "--primer_seq", required = False, help = "primer sequence to be used for read filtering")
    pre.add_argument("-barcode", "--barcode", required = False, help = "barcode table used for filtering on read2 (no header)")
    # alignment 
    align = sub_parsers.add_parser("align", help = "align fastq to reference genome", parents = [parser], add_help = False)
    align.add_argument("-fastq", "--fastq", nargs = 2, required = True, help="fastq input file (for paired-end reads, order in R1 R2)")
    align.add_argument("-ref", "--reference", required = False, help = "reference genome/genome index used for alignment")
    # screen for targets
    screen = sub_parsers.add_parser("screen", help = "screen alignment file", parents = [parser], add_help = False)
    screen.add_argument("-bam", "--bam", required = True, help = 'bam alignment file')
    screen.add_argument("-min_quality", "--min_quality", default = 10, type = int, help = "minimum mapping quality used for alignment file")
    screen.add_argument("-feature", "--feature", default = "CDS", help = "feature used to search for in GTF (default: CDS)")
    screen.add_argument("-frame", default = 0, type = int, help = "target position in frame", choices = [0, 1, 2])
    screen.add_argument("-donor", "--donor", required = False, help = "donor sequence used to filter read")
    screen.add_argument("-gRNA", "--gRNA", required = False, help = "guide RNA and its sequence feature")
    # parse arguments
    args=parser2.parse_args()
    return args

def parse_read1(read1_df:pd.DataFrame, primer_seq:str):
    """R1 is used to select reads containing specific SA (splicing acceptor)
    Return: a dataframe of read1 with designated primer [name, seq, code]
    """
    if primer_seq != "":
        read1_df["primer"] = read1_df["seq"].str.slice(0, len(primer_seq))
        read1_df_select = read1_df[read1_df["primer"] == primer_seq].drop(labels = "primer", axis = 1)
    else:
        read1_df_select = read1_df.copy()
    return read1_df_select

def parse_read2(read2_df:pd.DataFrame, barcode_list:list):
    """
    Return: a dataframe of read2 with with valid barcode [name, seq, code]
    """
    read2_df["barcode"] = read2_df["seq"].str.slice(0, len(barcode_list[0]))
    # the first 8bp is cell barcode (the barcode should match to cell barcode table given) 
    # it has the anticipated adapter and CCC cluster
    read2_df_select = read2_df[read2_df["barcode"].isin(barcode_list)].drop(labels = "barcode", axis = 1)
    return read2_df_select

def write_read(select_df:pd.DataFrame, output:str, start:int, end = None):
    select_df = select_df.copy()
    if start != 0:
        select_df["cut"] = select_df["seq"].str.slice(0, start)
        # update sequence for write into new file 
        select_df["seq"] = select_df["seq"].str.slice(start, end)
        select_df["quality"] = select_df["quality"].str.slice(start, end)
        with open(output, "w") as fr:
            for row in select_df.itertuples():
                fr.write(f"@{row.name} head:{row.cut} {row.comment}\n{row.seq}\n+\n{row.quality}\n")
    else:
        with open(output, "w") as fr:
            for row in select_df.itertuples():
                fr.write(f"@{row.name} {row.comment}\n{row.seq}\n+\n{row.quality}\n")

def STAR_align(r1, index_dir, gtf, output, n):
    command = f"STAR --runThreadN {n} --genomeDir {index_dir} --sjdbGTFfile {gtf} --readFilesIn {r1} --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix {output}/{output} --outSAMtype BAM SortedByCoordinate && samtools index -@ {n} {output}/{output}Aligned.sortedByCoord.out.bam"
    # disable twopassMode --twopassMode Basic 
    # --outFilterIntronMotifs RemoveNoncanonicalUnannotated (filter out alignments contain non-canonical unannotated junctions, annotated non-canonical junctions will be kept.)
    logging.info(f"Command\t{command}")
    subprocess.call(command, shell = True)
    # remove temporatory directories when exists
    import shutil 
    if os.path.exists(f"{output}/{output}_STARgenome"):
        try:
            shutil.rmtree(f"{output}/{output}_STARgenome")
        except Exception as e:
            print(e)
    if os.path.exists(f"{output}/{output}_STARpass1"):
        try:
            shutil.rmtree(f"{output}/{output}_STARpass1")
        except Exception as e:
            print(e)

def BWA_align(r1, ref, output, n):
    command = f"bwa mem -t {n} {ref} {r1} -L 0,0 | samtools view -Su - | samtools sort -@ {n} - -o {output}.bam && samtools index -@ {n} {output}.bam"
    logging.info(f"Command\t{command}")
    subprocess.call(command, shell = True)

def STAR_log(output):
    with open(os.path.join(output, output + "Log.final.out"), "r") as fr:
        for line in fr.readlines():
            line = line.strip()
            if "Number of input reads |" in line:
                input_r = line.split("Number of input reads |\t")[-1]
            if "Uniquely mapped reads number |" in line or "Uniquely mapped reads % |" in line:
                if "Uniquely mapped reads number |" in line:
                    unique_r = line.split("Uniquely mapped reads number |\t")[-1]
                if "Uniquely mapped reads % |" in line:
                    unique_p = line.split("Uniquely mapped reads % |\t")[-1].strip("%")
    logging.info(f"Align input\t{input_r}")
    logging.info(f"Align unique\t{unique_r}({unique_p}%)")

def BWA_log(align_read, output):
    fq_handle = pysam.FastxFile(align_read)
    total_read = 0
    for entry in fq_handle:
        total_read += 1
    logging.info(f"Align input\t{total_read}")
    read_bed = pd.read_table(os.path.join(output, output + ".bed"), sep = "\t", header = None)
    read_num = len(read_bed)
    read_p = round(read_num/total_read*100, 2)
    logging.info(f"Align unique\t{read_num}({read_p}%)")

def read_feature_fun(read_df, feature_df, label = "cDNA"):
    # find read overlap feature region of expectation 
    read_feature = bf.overlap(read_df, feature_df, how = "inner")
    # read no gap within 
    read_feature = read_feature[~read_feature["cigar"].str.contains("I|D")]
    read_feature = read_feature[read_feature["strand"] != read_feature["strand_"]]
    # check strandness of read and feature 
    if label.upper() == "CDNA":
        read_no_junction = read_feature[~(read_feature["cigar"].str.contains("N")) & (read_feature["start"] >= read_feature["start_"]) & (read_feature["end"] <= read_feature["end_"])]
        # across junction region
        read_other = read_feature[~read_feature["name"].isin(read_no_junction["name"])] 
        read_other.index = range(len(read_other))
        # check junction region (read junction match to exon junction)
        cigar_df = read_cigar(read_other)
        row_list = list()
        for row in read_other.itertuples():
            rowi = row.Index 
            r_cigar = list(zip(cigar_df.iloc[rowi]["char"], cigar_df.iloc[rowi]["digit"]))
            i = 0
            # when read across junction
            if "N" in row.cigar:
                for r in r_cigar:
                    i += int(r[-1])
                    if r[0] == "N":
                        if row.start + i == row.start_:
                            row_list.append(rowi)
                            break 
        read_junction = read_other.iloc[row_list]
        read_qualified = pd.concat([read_no_junction, read_junction], axis = 0)
        return read_qualified
    if label.upper() == "DNA":
        return read_feature

def read_frame_fun(read_qualified):
    """read should match exact one end of CDS and with appropriate orientation as opposing to gene strand"""
    read_qualified["junction_start"] = np.where((read_qualified["strand"] == "+") & (read_qualified["cigar_char"].str[0] == "S"), read_qualified["start"]-read_qualified["cigar_len"].str[0].astype(int), read_qualified["start"])
    read_qualified["junction_end"] = np.where((read_qualified["strand"] == "-") & (read_qualified["cigar_char"].str[-1] == "S"), read_qualified["end"]+read_qualified["cigar_len"].str[-1].astype(int), read_qualified["end"])
    #
    read_qualified["read_start"] = np.where(read_qualified["strand"] == "-", read_qualified["junction_end"], read_qualified["junction_start"])
    read_qualified["cds_1aa"] = np.where(read_qualified["strand_"] == "-", read_qualified["end_"]-read_qualified["frame_"], read_qualified["start_"]+read_qualified["frame_"])
    # 
    read_qualified["inframe"] = abs(read_qualified["cds_1aa"] - read_qualified["read_start"])%3
    #
    read_frame0 = read_qualified[read_qualified["inframe"] == 0].copy()
    read_frame1 = read_qualified[read_qualified["inframe"] == 3-1].copy()
    read_frame2 = read_qualified[read_qualified["inframe"] == 3-2].copy()
    return (read_frame0, read_frame1, read_frame2)

def frame_df_out(frame_df):
    # read_df
    # "chrom", "start", "end", "name", "score", "strand", "cigar"
    # feature_df
    # "chrom", "start", "end", "frame", "transcript", "gene", "symbol", "strand"
    frame_clean = frame_df.drop(["score", "chrom_", "start_", "end_", "frame_", "transcript_", "strand_", "char", "digit", "junction_start", "junction_end", "read_start", "cds_1aa", "inframe"], axis = 1)
    frame_clean = frame_clean.drop_duplicates(keep = "first")
    return frame_clean

def frame_count(frame_df):
    c_list = list()
    for g, gdf in frame_df.groupby(["gene_", "symbol_"]):
        c_list.append((g[0], g[1], len(set(gdf["name"])), list(set(gdf["chrom"]))[0], gdf["start"].min()))
    frame_out = pd.DataFrame(c_list, columns = ["gene", "symbol", "count", "chrom", "pos"])
    return frame_out

def read_donor_filter(read_featured, donor):
    from Bio import Align
    read_featured["clip_size"] = np.where((read_featured["strand"] == "+") & (read_featured["cigar_char"].str[0] == "S"), read_featured["cigar_len"].str[0].astype(int), np.where((read_featured["strand"] == "-") & (read_featured["cigar_char"].str[-1] == "S"), read_featured["cigar_len"].str[-1].astype(int), 0))
    # for read not reported as soft-clipping, I retrieve partial for alignment
    read_featured["clip_size"] = np.where(read_featured["clip_size"] < len(donor)*0.8, 16, read_featured["clip_size"])
    # 
    clip_seq_list = list()
    for row in read_featured.itertuples():
        if row.strand == "+":
            clip_seq = row.seq[:row.clip_size]
        else:
            clip_seq = reverse_complement(row.seq[-row.clip_size:])
        clip_seq_list.append(clip_seq)
    read_featured["clip_seq"] = clip_seq_list
    # 
    clip_df = read_featured[["clip_seq"]].drop_duplicates(keep = "first", ignore_index = True)
    # 
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -4
    aligner.internal_open_gap_score = -6
    aligner.internal_extend_gap_score = -1
    aligner.target_left_gap_score = -6
    aligner.target_right_gap_score = 0
    aligner.query_left_gap_score = -6
    aligner.query_right_gap_score = 0
    # 
    align_score = list()
    for row in clip_df.itertuples():
        for alignment in aligner.align(donor, row.clip_seq): # (target ,query)
            align_score.append(alignment.score)
            break
    clip_df["donor_score"] = align_score
    read_featured_scored = pd.merge(read_featured, clip_df, on = "clip_seq", how = "left")
    read_featured_scored = read_featured_scored[read_featured_scored["donor_score"] > 5].copy()
    return read_featured_scored


def main():
    args = args_parser()
    #print(args)
    try:
        output = args.output.strip("/")
        if not os.path.exists(output):
            print("Create directory ....")
            os.mkdir(output)
        if args.command == "pre":
            logging.basicConfig(filename = os.path.join(output, output+".log"), filemode = "w", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
        else:
            logging.basicConfig(filename = os.path.join(output, output+".log"), filemode = "a", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
        if args.assembly != None:
            logging.info(f"Assembly\t{args.assembly}")
    except Exception as e:
        print(e)
        exit(1)
    # requied argument
    # preprocessing enabled
    if args.command == "pre":
        r1, r2 = args.fastq
        r1_df = fq2df(r1)
        r2_df = fq2df(r2)
        logging.info(f"Total read\t{len(r1_df)}")
        logging.info(f"### Read cleaning ....")
        print("Read cleaning ....")
        primer_seq = args.primer_seq
        primer_seq = primer_seq.upper()
        logging.info(f"Primer seq\t{primer_seq}")
        #print(f"Primer: {primer_seq}")
        r1_head = len(primer_seq)
        r1_select_df = parse_read1(r1_df, primer_seq)
        if args.barcode != None:
            logging.info("Barcode mode\tTrue")
            barcode_df = pd.read_table(args.barcode, sep = "\t", header = None, names = ["barcode"])
            barcode_list = barcode_df["barcode"].tolist()
            r2_select_df = parse_read2(r2_df, barcode_list)
            r2_head = len(barcode_list[0])
            select_read = list(set(r1_select_df["name"]).intersection(set(r2_select_df["name"])))
        else:
            logging.info("Barcode mode\tFalse")
            r2_head = 0
            select_read = r1_select_df["name"].tolist()
        # filter read based on selection 
        logging.info(f"Clean Read\t{len(select_read)}({round(len(select_read)/len(r1_df)*100, 2)}%)")
        r1_select_df = r1_df[r1_df["name"].isin(select_read)]
        r2_select_df = r2_df[r2_df["name"].isin(select_read)]
        # write r1/r2 
        if args.read_length != None:
            r1_tail = r1_head + args.read_length
            r2_tail = r2_head + args.read_length
        else:
            r1_tail = None 
            r2_tail = None 
        print("Write new fastq ....")
        logging.info(f"sub read length\t{args.read_length}")
        write_read(r1_select_df, os.path.join(output, output + "_R1.fq"), start=r1_head, end=r1_tail)
        write_read(r2_select_df, os.path.join(output, output + "_R2.fq"), start=r2_head, end=r2_tail)
        align_r1 = os.path.join(output, output + "_R1.fq")
        align_r2 = os.path.join(output, output + "_R2.fq")
    # when alignment enabled 
    if args.command == "align":
        logging.info(f"Seq data\t{args.dtype}")
        try:
            r1, r2 = args.fastq
            align_r1 = r1
            align_r2 = r2
            align_index = args.reference
            if args.dtype.upper() == "CDNA":
                gtf = args.gtf
                if gtf.endswith(".gz"):
                    gtf = gtf.split(".gz")[0]
                gtf_gz = gtf + ".gz"
                if not os.path.exists(os.path.join(output, output + "Log.final.out")):
                    logging.info("### Read alignment ....")
                    print("Read alignment ....")
                    STAR_align(align_r1, align_index, gtf, output, n = 4)
                bam2bed(os.path.join(output, output + "Aligned.sortedByCoord.out.bam"), label = "CDNA").to_csv(os.path.join(output, output+".bed"), sep = "\t", header = False, index = False)
                print("Alignment report ....")
                logging.info("### Alignment report ....")
                STAR_log(output)
            if args.dtype.upper() == "DNA":
                if not os.path.exists(os.path.join(output, output + ".bam")):
                    logging.info("### Read alignment ....")
                    BWA_align(align_r1, align_index, os.path.join(output, output), n = 4)
                bam2bed(os.path.join(output, output + ".bam"), label = "DNA").to_csv(os.path.join(output, output+".bed"), sep = "\t", header = False, index = False)
                print("Alignment reprot .....")
                logging.info("### Alignment report ....")
                BWA_log(align_r1, output)
            if args.dtype.upper() not in ["CDNA", "DNA"]:
                print("Sequencing datatype is not recognized.")
                exit(1)
        except Exception as e:
            print(e)
    # when screen enabled 
    if args.command == "screen":
        logging.info("### Parse alignment ....")
        if args.gtf == None:
            print("Please provide gtf file")
            exit(1)
        gtf = args.gtf
        if gtf.endswith(".gz"):
            gtf = gtf.split(".gz")[0]
        gtf_gz = gtf + ".gz"
        print("Read gtf feature ....")
        feature_df = gtf_parse(gtf_gz, args.feature)
        # read bed file 
        read_bed = pd.read_table(os.path.join(output, output + ".bed"), sep = "\t", header = None, names = ["chrom", "start", "end", "name", "score", "strand", "cigar", "seq"])
        # add cigar_char, cigar_len to read bed file
        read_bed = pd.merge(read_bed, cigar_soft(read_bed), left_index = True, right_index = True)
        ######### filtering bed by mapping quality 
        read_clean = read_bed[read_bed["score"] > args.min_quality].copy()
        # total number of reads (only for reads with unique/primary alignment, and score > cutoff)
        bed_num = len(set(read_clean["name"]))
        logging.info(f"Total alignment\t{bed_num}")
        # filter bed by overlapping feature
        print("Filter feature ....")
        read_featured = read_feature_fun(read_clean, feature_df, label = args.dtype)
        read_featured_num = len(set(read_featured["name"]))
        logging.info(f"{args.feature} alignment\t{read_featured_num}({round(read_featured_num/bed_num*100, 2)}%)")
        genes_featured_num = len(set(read_featured["gene_"]))
        logging.info(f"Genes w/ read\t{genes_featured_num}")
        # filter bed by donor sequence
        if args.donor != None:
            logging.info(f"Donor seq to check\t{args.donor}")
            read_donored = read_donor_filter(read_featured, args.donor)
            read_donored_num = len(set(read_donored["name"]))
            logging.info(f"donor filter\t{read_donored_num}({round(read_donored_num/bed_num*100, 2)}%)")
            genes_donored_num = len(set(read_donored["gene_"]))
            logging.info(f"Genes w/ donor-containing-read\t{genes_donored_num}")
        else:
            read_donored = read_featured
        # frame filtering
        print("Frame filtering ....")
        logging.info(f"Frame\t{args.frame}")
        frame0_df, frame1_df, frame2_df = read_frame_fun(read_donored)
        frame0_gene_num = len(set(frame0_df["gene_"]))
        frame1_gene_num = len(set(frame1_df["gene_"]))
        frame2_gene_num = len(set(frame2_df["gene_"]))
        logging.info(f"Frame0 gene#\t{frame0_gene_num}")
        logging.info(f"Frame1 gene#\t{frame1_gene_num}")
        logging.info(f"Frame2 gene#\t{frame2_gene_num}")
        if args.frame == 0:
            bamfilter(args.bam, list(set(frame0_df["name"])), output)
            target_df = frame_count(frame0_df)
            target_df.to_csv(os.path.join(output, output + "_inframe.txt"), sep = "\t", header = True, index = False)
        elif args.frame == 1:
            bamfilter(args.bam, list(set(frame1_df["name"])), output)
            target_df = frame_count(frame1_df)
            target_df.to_csv(os.path.join(output, output + "_inframe.txt"), sep = "\t", header = True, index = False)
        elif args.frame == 2:
            bamfilter(args.bam, list(set(frame2_df["name"])), output)
            target_df = frame_count(frame2_df)
            target_df.to_csv(os.path.join(output, output + "_inframe.txt"), sep = "\t", header = True, index = False)
        else:
            print("Frame is not recognized")
            exit(1)
        # filter for genes of target
        # get target candidates
        if args.gRNA != None:
            print("Read target candidates ....")
            gRNA_gene_df = pd.read_table(args.gRNA, sep = "\t", header = 0)
            target_num = len(set(gRNA_gene_df["symbol"]))
            logging.info(f"Target candidate\t{target_num}")
        else:
            exit(0)
        expected_gene = set(gRNA_gene_df["symbol"])
        observed_gene = set(target_df["symbol"])
        on_target = observed_gene.intersection(expected_gene)
        out_of_target = observed_gene.difference(expected_gene)
        logging.info(f"On-target\t{len(on_target)}")
        logging.info(f"Out-of-target\t{len(out_of_target)}")

if __name__ == "__main__":
    main()
