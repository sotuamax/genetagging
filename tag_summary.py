#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os 
import sys
import pandas as pd
import numpy as np
import argparse


"""
summarize GFP results 
filtering on GFP results based on read-counts 

"""
def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter, description="Use bed file output from tag_align.py (overlap between read and gene annotation feature), to report genes tagged by GFP")
    # important parameters 
    parser.add_argument("-sample", "--sample", help="sample prefix used to find its bed file")
    parser.add_argument("-outdir", "--outdir", default = ".", help="output directory name")
    parser.add_argument("-bc", "--bc", action = "store_true", help = "enable barcode mode. ")
    # parser.add_argument("-p", "--plot", help = "plot GFP read distribution", action = "store_true")
    parser.add_argument("-frame", default = 0, type = int, help = "target position in frame", choices = [0, 1, 2])
    args=parser.parse_args()
    return args

def main(): 
    args = args_parser()
    sample = args.sample
    output_dir = args.outdir
    frame = args.frame
    # get preprocessing info 
    # try:
    #     print("Load preprocessing data ... ")
    #     # pre_log = pd.read_table(os.path.join(output_dir, sample + ".pre.txt"), sep = "\t", header = 0)
    #     # raw_qc = pd.read_table(os.path.join(output_dir, sample + ".stat"), sep = "\t", header = 0)
    # except Exception as e:
    #     print(e)
    #     exit(1)
    # get processed bed for good quality reads on CDS (NO frame filtering yet)
    print("Load read alignment bed ...")
    read_feature = pd.read_table(os.path.join(output_dir, sample + ".qc.bed"), sep = "\t", header = 0)
    if len(set(read_feature["name"])) != len(read_feature):
        raise ValueError("Read name are not unqiue in the bed file.")
    print("Read frame stats ...")
    from utilities.excel_tools import excel_write
    if args.bc:
        read_feature["barcode"] = read_feature["name"].str.split(":", expand = True).iloc[:, -1]
        read_feature["index"] = read_feature["barcode"].str.slice(0,15)
        # summarize library read count
        index_CDS = read_feature.groupby("index").size().reset_index(name = "CDS_align")
        index_frame = read_feature.groupby(["index", "frame_"]).size().reset_index(name = "frame_count")
        index_CDS_frame = pd.merge(index_CDS, index_frame, on = "index")
        index_CDS_frame["frame_yield"] = round(index_CDS_frame["frame_count"]/index_CDS_frame["CDS_align"]*100)
        index_CDS_frame.to_csv(os.path.join(output_dir, sample + ".frame.stat"), sep = "\t", header = True, index = False)
        # combine preprocessing step with CDS align/frame yield 
        print("Combine preprocessing log with frame info ...")
        #pre_combine = pd.merge(pre_log, index_CDS_frame.query("frame_ == @frame"), on = "index")
        #pre_combine.sort_values("library", inplace = True, ignore_index=True)
        #excel_write(pre_combine, os.path.join(output_dir, sample + ".pre.xlsx"))
        # frame_stat = pd.merge(CDS_read, frame_read, on = "index")
        # frame_stat["frame_ratio"] = round(frame_stat['frame_read']/frame_stat['CDS_read']*100)
        # frame_stat.to_csv(os.path.join(output_dir, sample + ".frame"), sep = "\t", header = True, index = False)
        print("Collect cell level frame info ...")
        barcode_read = read_feature.groupby(["barcode"]).size().reset_index(name = "barcode_read")
        frame_read = read_feature.groupby(["barcode", "frame_"]).size().reset_index(name = "frame_read")
        barcode_frame_stat = pd.merge(barcode_read, frame_read, on = "barcode")
        barcode_frame_stat["frame_ratio"] = round(barcode_frame_stat["frame_read"]/barcode_frame_stat["barcode_read"]*100)
        gene_read = read_feature.groupby(["barcode", "frame_", "gene_", "symbol_"]).size().reset_index(name = "gene_frame_read")
        barcode_frame_stat = pd.merge(barcode_frame_stat, gene_read, on = ["barcode", "frame_"])
        barcode_frame_stat["gene_frame_ratio"] = round(barcode_frame_stat["gene_frame_read"]/barcode_frame_stat["frame_read"]*100)
        barcode_frame_stat.to_csv(os.path.join(output_dir, sample + ".gene.stat"), sep = "\t", header = True, index = False)
        gene_clean = barcode_frame_stat.query("gene_frame_read >= 3 and gene_frame_ratio > 80 and frame_ratio > 80")
        gene_clean.to_csv(os.path.join(output_dir, sample + ".gene.clean.stat"), sep = "\t", header = True, index = False)
        target_frame_gene = gene_clean.query("frame_ == @frame")
        print(f"Identified genes with target frame {frame}: {len(target_frame_gene)}")
        print(f"target frame {frame} ratio: {round(len(target_frame_gene)/len(gene_clean)*100)}%")
        print("Report genes ...")
        excel_write(target_frame_gene, os.path.join(output_dir, sample + ".gene.xlsx"))
    else:
        excel_write(read_feature, os.path.join(output_dir, sample + ".pre.xlsx"))
        read_inframe = read_feature.query("frame_ == @frame").copy()
        excel_write(read_inframe, os.path.join(output_dir, sample + ".inframe.xlsx"))
        gene_uni = len(set(read_inframe["gene_"]))
        print(f"Identified genes with target frame {frame}: {gene_uni} ")
        print(f"target frame {frame} ratio: {round(len(read_inframe)/len(read_feature)*100, 2)}%")
        print("Report genes ... ")
        read_inframe["fusion"] = np.where(read_inframe["strand"] == "+", read_inframe["start"], read_inframe["end"])
        target_frame_gene = read_inframe.groupby("gene_").size().reset_index(name = "read_count")
        gene_fusion_site = read_inframe.groupby("gene_")["fusion"].nunique().reset_index(name = "fusion_number")
        target_gene = pd.merge(pd.merge(target_frame_gene, gene_fusion_site, on = "gene_"), read_inframe[["gene_", "symbol_"]].drop_duplicates(keep = "first"), on = "gene_")
        #target_frame_gene[target_frame_gene["read_count"] >= 3]
        #excel_write(target_gene, os.path.join(output_dir, sample + ".gene.xlsx"))
        target_gene.to_csv(os.path.join(output_dir, sample + ".gene.csv"), sep = "\t", header = True, index = False)
        fusion, fusion_freq =np.unique(target_gene["fusion_number"], return_counts=True)
        fusion_freq_df = pd.DataFrame({"fusion_site":fusion, "freq":fusion_freq})
        fusion_freq_df["prop"] = fusion_freq_df["freq"]/sum(fusion_freq_df["freq"])*100
        fusion_freq_df.to_csv(os.path.join(output_dir, sample + ".fusion.csv"), sep = "\t", header = True, index = False)
        exit(0)
    # re-write results into a master excel 
    # pre_excel = pd.read_excel(os.path.join(output_dir, sample + ".pre.xlsx"))
    gene_excel = pd.read_excel(os.path.join(output_dir, sample + ".gene.xlsx"))
    gene_excel["index"] = gene_excel["barcode"].str.slice(0,15)
    # gene_excel = pd.merge(pre_excel[["library", "index"]], gene_excel, how = "right", on = "index")
    # i_list = [[i, len(i_df["gene_"]), len(set(i_df["gene_"]))] for i,i_df in gene_excel.groupby("index")]
    # gene_count = pd.DataFrame(i_list, columns = ["index", "gene", "unique_gene"])
    # pre_excel = pd.merge(pre_excel, gene_count, on = "index")
    writer = pd.ExcelWriter(os.path.join(output_dir, sample + ".xlsx"), mode = "w")
    # pre_excel.to_excel(writer, sheet_name = "summary", index = False)
    gene_excel.to_excel(writer, sheet_name = "flag_gene", index = False)
    writer.close()


if __name__ == "__main__":
    main()
