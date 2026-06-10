#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program is to clean sequencing data based on the primer given (R1, enriched), 
and assign cell identify based on cell barcode given in R2.
Output file is clean and cell barcode labeled R1 fastq file 
1). start with primer; 
2). valid cell barcode. 

tag_align.py -sample GC760804GFP -seq_type cDNA -outdir output/GC760804 -n 24 -ref $hg38_star -qc -report -gtf $gtf 

For ChiC data: 
tag_align.py -sample $sample -seq_type DNA -outdir $output_dir -n 32 -ref $hg38 -report -qc 
For GFP data:
tag_align.py -sample ${sample}GFP -seq_type cDNA -outdir $output_dir -n 32 -ref $hg38_star -qc -gtf $gtf_pick

"""

import os 
import sys
import pandas as pd
import bioframe as bf 
import numpy as np
import argparse
import pysam 
import subprocess
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
from utilities.bam_tools import flagstats
from utilities.bam_tools import bamfilter 
from utilities.bam_tools import bam2bw
import matplotlib.pyplot as plt

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("-sample", "--sample", help="sample prefix used to find fastq file")
    parser.add_argument("-outdir", "--outdir", default = ".", help="output directory name (also for sample fastq location)" \
                        "output files: " \
                        "*_barcode.stat -> ChiC read count at single cell level; " \
                        "*_clean_barcode.stat -> ChiC read count filtered by GFP at single cell level;" \
                        "*GFP.stat -> GFP barcode statistics; "
                        "*GFP.gene.stat -> GFP genes unfiltered; " \
                        "*GFP.gene.clean.stat -> GFP genes filtered by frame and single target (dominate by one gene)")
    parser.add_argument("-n", "--threads", help = "number of threads for alignment", default = 1, type = int)
    parser.add_argument("-seq_type", "--seq_type", help = "sequencing data type", choices = ["cDNA", "DNA", "merge"])
    parser.add_argument("-ref", "--reference", required=False, help = "reference genome index used for alignment (note that for cDNA, star index required)")
    parser.add_argument("-dedup", "--dedup", action = "store_true", help = "deduplicate on ChiC alignment")
    parser.add_argument("-qc", action = "store_true", help = "perform QC on BAM alignment file (for DNA/cDNA)")
    parser.add_argument("-report", action = "store_true", help = "generate DNA report for read count on cell barcode, "
    "or generate cDNA report for cDNA and its overlap with CDS site (only need to run it once)")
    parser.add_argument("-fusion", "--fusion", action = "store_true", help = "identify GFP fusion genes (only need to run it once)")
    parser.add_argument("-frame", default = 0, type = int, help = "target position in frame", choices = [0, 1, 2])
    parser.add_argument("-f", "--force", action = "store_true", help = "force alignment even when bam file exists")
    parser.add_argument("-gtf", "--gtf", required = False, help = "reference gtf file used for RNA-seq alignment (for cDNA)")
    parser.add_argument("-GFP_filter", "--GFP_filter", action = "store_true", help = "Filter ChiC data based on GFP barcode (must after GFP report). ")
    parser.add_argument("-cap", "--cap", help = "CAP file", required = False)
    args=parser.parse_args()
    return args

def mode_freq(x):
    return x.value_counts().iloc[0]

def main():
    args = args_parser()
    sample = args.sample
    output_dir = args.outdir
    bam_file = os.path.join(output_dir, sample + ".bam")
    n= args.threads
    try:
        os.mkdir(output_dir)
        print(f"Create {output_dir}")
    except:
        pass
    ### align for DNA seq
    if args.seq_type == "DNA":
        from utilities.align_tools import BWA_SE
        if not os.path.exists(bam_file) or args.force:
            align_r1 = os.path.join(output_dir, sample + "_R1.fastq.gz")
            if not os.path.exists(align_r1):
                raise IOError("Cannot find R1 for alignment.")
            print(f"Locate {sample} file ...")
            print(f"Perform {args.seq_type} alignment ...")
            # perform single-end BWA alignment 
            BWA_SE(align_r1, args.reference, os.path.join(output_dir, args.sample), n)
            flagstats(bam_file, n)
        if args.qc:
            # check soft-clip side at 5' end for QC 
            qc_bam = bam_file.replace(".bam", ".qc.bam")
            #if not os.path.exists(qc_bam):
            print(f"Process {bam_file} ...")
            # bamfilter in default, only retain reads on chromosome region (no ChrM/chrY)
            # perform QC on BAM file 
            if not os.path.exists(qc_bam):
                print("Perform QC ...")
                bamfilter(bam_file, qc_bam, clip_check=True, threads = n, chrom=True)
                flagstats(qc_bam, n)
        if args.dedup:
            # enable single-cell deduplicate process
            from utilities.bam_tools import sc_markdup
            print(f"Process {bam_file} ... ")
            # markdup(bam_file, os.path.join(output_dir, args.sample) + ".dedup")
            # make sure that single cell deduplication is on cell barcode and coordinates and strand
            if os.path.exists(bam_file.replace(".bam", ".qc.bam")):
                bam_file = bam_file.replace(".bam", ".qc.bam")
            dedup_bam = bam_file.replace(".bam", ".dedup.bam")
            bw = dedup_bam.replace(".bam", ".bw")
            if not os.path.exists(dedup_bam):
                print("Deduplicate on single cell level ...")
                sc_markdup(bam_file, dedup_bam, n)
                flagstats(dedup_bam, n)
            # if not os.path.exists(bw):
            #     bam2bw(qc_bam, bw, n)
        if args.report:
            bam_file = os.path.join(output_dir, sample + ".dedup.qc.bam")
            if os.path.exists(bam_file):
                pass 
            else:
                bam_file = os.path.join(output_dir, sample + ".qc.dedup.bam")
            bam_handle = pysam.AlignmentFile(bam_file, "rb", threads = n)
            if not os.path.exists(os.path.join(output_dir, sample + "_dna.stat")):
                print(f"Generate report for cell barcode on {sample}...")
                print(f"Process {os.path.basename(bam_file)} ... ")
                read_barcode_count = dict()
                for read in bam_handle.fetch():
                    if (read.query_name).split(":")[-1] not in read_barcode_count:
                        read_barcode_count[(read.query_name).split(":")[-1]] = 0
                    read_barcode_count[(read.query_name).split(":")[-1]] += 1
                #read_barcode_list = [ for read in bam_handle.fetch()]
                #barcode, barcode_freq = np.unique(read_barcode_list, return_counts = True)
                barcode_freq_df = pd.DataFrame.from_dict(read_barcode_count, orient="index").reset_index()
                barcode_freq_df.columns = ["barcode", "freq"]
                barcode_freq_df.to_csv(os.path.join(output_dir, sample + "_dna.stat"), sep = "\t", header = True, index = False)
            else:
                barcode_freq_df = pd.read_table(os.path.join(output_dir, sample + "_dna.stat"), sep = "\t", header = 0)
            barcode_freq_df.sort_values(by = "freq", ascending=False, inplace = True, ignore_index = True)
            # make plots for ChiC read per cell
            fig, axes = plt.subplots(1, 2, figsize = (10,4))
            axes[0].plot(range(len(barcode_freq_df)), barcode_freq_df["freq"], "o", markersize = 0.2)
            axes[0].set_xscale('log')
            axes[0].set_yscale('log')
            axes[0].set_xlabel("cell barcode ranked by read count")
            axes[0].set_ylabel("read count per cell barcode")
            axes[1].hist(np.log10(barcode_freq_df["freq"]), bins = 50, orientation='horizontal')
            axes[1].set_ylabel("log(read count per cell barcode)")
            axes[1].set_xlabel("freq")
            plt.savefig(os.path.join(output_dir, sample + "_dna.png"), dpi = 300)
            # combine with GFP (when GFP data available)
            # ChiC read per cell for cells with valid GFP
            # if os.path.exists(os.path.join(output_dir, sample + "GFP.xlsx")):
            #     print("Generate report for cell barcode w/ GFP filtering ...")
            #     tag_df = pd.read_excel(os.path.join(output_dir, sample + "GFP.xlsx"), sheet_name  = "flag_gene")
            #     if "library" not in tag_df.columns:
            #         library_df = pd.read_excel(os.path.join(output_dir, sample + "GFP.xlsx"), sheet_name  = "summary")
            #         tag_df["index"] = tag_df["barcode"].str.slice(0,15)
            #         tag_df = pd.merge(library_df[["library", "index"]].drop_duplicates(keep = "first"), tag_df, on = "index", how = "right")
            #     tag_df.columns = [c.strip("_") for c in tag_df.columns]
            #     # ChiC_barcode_freq = pd.read_table(os.path.join(output_dir, sample + "_barcode.stat"), sep = "\t", header = 0)
            #     # ChiC_barcode_freq.sort_values(by = ['freq'], ascending=False, inplace = True, ignore_index = True)
            #     gfp_tag_copy = tag_df[["library", "index", "barcode", "gene", "symbol", "gene_frame_read"]].copy()
            #     gfp_tag_copy["i7"] = gfp_tag_copy['barcode'].str.split("+", expand = True)[0]
            #     barcode_freq_df["i7"] = barcode_freq_df["barcode"].str.split("+", expand = True)[0]
            #     gfp_tag_copy["b"] = gfp_tag_copy["barcode"].str.split("+", expand = True)[2]
            #     barcode_freq_df["b"] = barcode_freq_df["barcode"].str.split("+", expand = True)[2]
            #     GFP_ChiC_df = pd.merge(gfp_tag_copy.rename(columns = {"barcode":"GFP_barcode"}), barcode_freq_df.rename(columns = {"barcode":"ChiC_barcode"}), on = ["b", "i7"], how = "inner")
            #     GFP_ChiC_df.drop(["i7", "b"], axis = 1, inplace = True)
            #     GFP_ChiC_df.rename(columns = {"gene_frame_read":"GFP_read_count", "freq":"ChiC_read_count"}, inplace = True)
            #     GFP_ChiC_df.sort_values(["ChiC_read_count"], ascending=[False], inplace = True, ignore_index=True)
            #     GFP_ChiC_df.to_csv(os.path.join(output_dir, sample + "_clean_barcode.stat"), sep = "\t", header = True, index = False)
            #     if args.GFP_filter:
            #         print("Filter ChiC BAM on GFP selected barcodes ...")
            #         GFP_ChiC_df = pd.read_table(os.path.join(output_dir, sample + "_clean_barcode.stat"), sep = "\t", header = 0)
            #         ChiC_barcode = GFP_ChiC_df["ChiC_barcode"].to_list()
            #         new_bam = bam_file.replace(".bam", ".GFP_filter.bam")
            #         if not os.path.exists(new_bam):
            #             with pysam.AlignmentFile(new_bam, "wb", template = bam_handle, threads = n) as newbam:
            #                 for read in bam_handle.fetch():
            #                     if read.query_name.split(":")[-1] in ChiC_barcode:
            #                         newbam.write(read)
            #         subprocess.call(f"samtools index -@ {n} {new_bam}", shell = True)
            #         new_bw = new_bam.replace(".bam", ".bw")
            #         subprocess.call(f"bamCoverage -b {new_bam} -o {new_bw} -of bigwig -p {n} --normalizeUsing None", shell = True)
                # enrich 
                # enrich_df = pd.read_table(os.path.join(output_dir, sample + ".enrich.txt"), sep = "\t", header = 0)
                # enrich_sum_df = pd.merge(enrich_df.groupby("barcode")["total"].sum().reset_index(), enrich_df.groupby("barcode")["overlap_peak"].sum().reset_index(), on ="barcode")
                # chrM_count = enrich_df.query("chrom == 'chrM'")[["barcode", "total"]].rename(columns = {"total":"chrM"})
                # enrich_sum_df = pd.merge(enrich_sum_df, chrM_count, on = "barcode")
                # enrich_sum_df.rename(columns = {"barcode":"ChiC_barcode"}, inplace = True)
                # enrich_sum_df["enrichment"] = enrich_sum_df["overlap_peak"]/enrich_sum_df["total"]
                # enrich_sum_df["chrM_ratio"] = enrich_sum_df["chrM"]/enrich_sum_df["total"]
                # GFP_ChiC_df = pd.merge(GFP_ChiC_df, enrich_sum_df, on = "ChiC_barcode")
                # GFP_ChiC_df["ChiC_read_count_log"] = np.log10(GFP_ChiC_df["ChiC_read_count"])
                # GFP_ChiC_df["GFP_read_count_log"] = np.log10(GFP_ChiC_df["GFP_read_count"])
                # GFP_ChiC = GFP_ChiC_df[["ChiC_read_count_log", "GFP_read_count_log", "enrichment", "chrM_ratio"]]
                # import seaborn as sns
                # pair = sns.pairplot(GFP_ChiC, plot_kws={'s': 5}, markers = ".")
                # pair.figure.savefig(os.path.join(output_dir, sample + "pair.png"), dpi=300, bbox_inches='tight')
                ### using GM (Gaussian Mixture model to separate cells and find the real single cell)
                # plot the relationship between GFP and ChiC read count 
                # from sklearn.mixture import GaussianMixture
                # gmm = GaussianMixture(n_components=2, random_state=0)
                # data = np.array(GFP_ChiC_df["ChiC_read_count_log"]).reshape(-1,1)
                # gmm.fit(data)
                # # gmm.predict(data)
                # # i = np.argmax(gmm.means_) # get the index for max-mean 
                # GFP_ChiC_df["label"] = -1
                # gmm_parameters = sorted([(gmm.means_[i][0],np.sqrt(gmm.covariances_[i][0][0])) for i in range(2)])
                # # default label == -1, group with smaller mean label == 0, group with larger mean label == 1 (1 is what we can use for further analysis)
                # for index, (m,sd) in enumerate(gmm_parameters):
                #     min_i = m-sd*1.645; max_i = m+sd*1.645
                #     GFP_ChiC_df["label"] = np.where((GFP_ChiC_df["ChiC_read_count_log"] > min_i) & (GFP_ChiC_df["ChiC_read_count_log"] < max_i), index, GFP_ChiC_df["label"])
                # # GFP_ChiC_df["label"] = np.where(GFP_ChiC_df["chrM_ratio"] < 0.2, GFP_ChiC_df["label"], -1)
                # ### make plots of the selected cell
                # fig, axes = plt.subplots(1, 2, figsize = (10, 4))
                # color_map = {0: 'red', 1:"red", -1: 'grey'}
                # colors = GFP_ChiC_df['label'].map(color_map)
                # axes[0].scatter(range(len(GFP_ChiC_df)), GFP_ChiC_df["ChiC_read_count"], c=colors, s = 0.2)
                # # .plot(, "o", markersize = 0.2, )
                # axes[0].set_xscale('log')
                # axes[0].set_yscale('log')
                # axes[0].set_xlabel("cell ranked by read count")
                # axes[0].set_ylabel("ChiC read count per cell")
                # axes[1].hist(GFP_ChiC_df["ChiC_read_count_log"], bins = 50, orientation='horizontal')
                # for i in range(2):
                #     m = gmm.means_[i][0]
                #     sd = np.sqrt(gmm.covariances_[i][0][0])
                #     min_i = m-sd*1.645; max_i = m+sd*1.645
                #     axes[1].axhline(y=min_i, color='black', linestyle='--', linewidth=1)
                #     axes[1].axhline(y=max_i, color='black', linestyle='--', linewidth=1)
                # axes[1].set_ylabel("log10(ChiC read count per cell)")
                # axes[1].set_xlabel("freq")
                # plt.savefig(os.path.join(output_dir, sample + "_clean_barcode.png"), dpi = 300)
                # GFP_ChiC_df.drop(["index", "ChiC_read_count_log", "GFP_read_count_log"], axis = 1).to_csv(os.path.join(output_dir, sample + "_clean_barcode.stat"), sep = "\t", header = True, index = False)
    ### align for RNA seq 
    if args.seq_type == "cDNA":
        from utilities.align_tools import STAR_SE
        gtf = args.gtf
        if gtf.endswith(".gz"):
            gtf = gtf.split(".gz")[0]
        ### 
        if not os.path.exists(bam_file) or args.force:
            align_r1 = os.path.join(output_dir, sample + "_R1.fastq.gz")
            if not os.path.exists(align_r1):
                raise IOError("Cannot find R1 for alignment.")
            print(f"Locate {sample} file ...")
            print(f"Perform {args.seq_type} alignment ...")
            STAR_SE(align_r1, args.reference, gtf, os.path.join(output_dir, sample), n = n, enable_novel = False)
            flagstats(bam_file, n)
        bw = bam_file.replace(".bam", ".bw")
        # if not os.path.exists(bw):
        #     bam2bw(bam_file, bw, n)
        if args.qc:
            # print("Perform QC on GFP BAM ...")
            # check soft-clip side at 5' end for QC 
            # from utilities.bam_tools import bamfilter 
            bam_handle = pysam.AlignmentFile(bam_file, "rb", threads = n)
            qc_bam = bam_file.replace(".bam", ".qc.bam")
            import re 
            if not os.path.exists(qc_bam):
                print("Perform QC ...")
                print(f"Process {bam_file} ...")
                with pysam.AlignmentFile(qc_bam, "wb", template = bam_handle, threads = n) as newbam:
                    for read in bam_handle.fetch():
                        # uniquely aligned and no mismatches 
                        if read.get_tag("NH") == 1 and read.get_tag("nM") == 0:
                            if "S" in read.cigarstring:
                                if (re.findall(r'[A-Z=]', read.cigarstring)[0] == "S" and read.is_forward) or (re.findall(r'[A-Z=]', read.cigarstring)[-1] == "S" and read.is_reverse) or ("I" in read.cigarstring) or ("D" in read.cigarstring):
                                    pass
                                else:
                                    newbam.write(read)
                            else:
                                newbam.write(read)
                subprocess.call(f"samtools index -@ {n} {qc_bam}", shell = True)
                flagstats(qc_bam, n)
            bw = qc_bam.replace(".bam", ".bw")
            # if not os.path.exists(bw):
            #     bam2bw(qc_bam, bw, n)
            qc_bam_handle = pysam.AlignmentFile(qc_bam, "rb", threads = n)
            barcode_dict = dict()
            for read in qc_bam_handle.fetch():
                if read.query_name.split(":")[-1] not in barcode_dict:
                    barcode_dict[read.query_name.split(":")[-1]] = 0
                barcode_dict[read.query_name.split(":")[-1]] += 1
            barcode_count_df = pd.DataFrame.from_dict(barcode_dict, orient = "index", columns = ["total_GFP"]).reset_index().rename(columns = {"index":"barcode"})
            barcode_count_df.to_csv(os.path.join(output_dir, sample + "_rna.stat"), sep = "\t", header = True, index = False)
            print("generate plot ...")
            barcode_count_df.sort_values(by = "total_GFP", ascending=False, inplace = True, ignore_index=True)
            fig, axes = plt.subplots(1, 2, figsize = (10,4))
            axes[0].plot(range(len(barcode_count_df)), barcode_count_df["total_GFP"], "o", markersize = 0.2)
            axes[0].set_xscale('log')
            axes[0].set_yscale('log')
            axes[0].set_xlabel("cell barcode ranked by GFP read count")
            axes[0].set_ylabel("GFP read count per cell barcode")
            axes[1].hist(np.log10(barcode_count_df["total_GFP"]), bins = 50, orientation='horizontal', color = "forestgreen")
            axes[1].set_ylabel("log(GFP read count per cell barcode)")
            axes[1].set_xlabel("freq")
            plt.savefig(os.path.join(output_dir, sample + "GFP_rna.png"), dpi = 300)
        if args.report: # report the overlap between read and CDS features in bed format 
            print("report the overlap between read and CDS ...")
            if os.path.exists(os.path.join(output_dir, sample + ".qc.bam")):
                bam_file = os.path.join(output_dir, sample + ".qc.bam")
            from utilities.gtf_tools import parse_gtf
            gtf_gz = gtf + ".gz"
            print("Parse CDS region ...")
            feature_df = parse_gtf(gtf_gz, "CDS").drop("transcript", axis = 1).drop_duplicates(keep = "first", ignore_index=True)
            # annotated frame is for 5' end 
            # update frame on 3' end (where donor seq to add)
            feature_df["frame"] = (feature_df["end"] - feature_df["frame"] - feature_df["start"])%3
            from utilities.bam_tools import bam2bed
            #filtered_list = list()
            print("Join alignment and CDS ...")
            #tmp_bam = bam_file.replace(".bam", ".tmp.bam")
            tmp_bed = bam_file.replace(".bam", ".bed")
            # delete all previous writing if exists
            # write chrom by chrom
            i = 0
            for c,feature_c in feature_df.groupby("chrom"):
                c_read = bam2bed(bam_file, chrom = c)
                if c_read is not None:
                    read_feature_overlap = bf.overlap(c_read, feature_c, how = "inner")
                    # read_feature_overlap.drop_duplicates(keep = "first", inplace = True, ignore_index = True)
                    # read_feature_overlap = read_feature_overlap[~read_feature_overlap["cigar"].str.contains("N")]
                    ### orientation & end overlap
                    scenario1 = (read_feature_overlap["strand"] == "+") & (read_feature_overlap['start'] == read_feature_overlap["start_"]) & (read_feature_overlap["strand_"] == "-")
                    scenario2 = (read_feature_overlap["strand"] == "-") & (read_feature_overlap['end'] == read_feature_overlap["end_"]) & (read_feature_overlap["strand_"] == "+")
                    read_feature_overlap_filtered = read_feature_overlap[scenario1 | scenario2]
                    ### per read match to single frame 
                    single_frame = read_feature_overlap_filtered.groupby("name")["frame_"].transform("nunique") == 1
                    read_feature_overlap_filtered = read_feature_overlap_filtered[single_frame]
                    ### per read (labeled by its name) match to single target gene
                    single_target = read_feature_overlap_filtered.groupby("name")["gene_"].transform("nunique") == 1
                    read_feature_overlap_filtered = read_feature_overlap_filtered[single_target].copy()
                    read_feature_overlap_filtered.drop(["chrom_", "start_", "end_", "strand_"], axis = 1, inplace = True)
                    read_feature_overlap_filtered.drop_duplicates(keep = "first", ignore_index = True, inplace = True)
                    if i == 0:
                        read_feature_overlap_filtered.to_csv(tmp_bed, mode = "w", sep = "\t", header = True, index = False)
                    else:
                        read_feature_overlap_filtered.to_csv(tmp_bed, mode = "a", sep = "\t", header = False, index = False)
                    i += 1
            print("Finish writing into bed file.")
        if args.fusion:
            print(f"Identify fusion genes for {sample} ...") # based on report output bed file
            # frame = args.frame
            if os.path.exists(os.path.join(output_dir, sample + ".qc.bam")):
                bam_file = os.path.join(output_dir, sample + ".qc.bam")
            bed = bam_file.replace(".bam", ".bed")
            if not os.path.exists(bed):
                print("No bed file identified")
                exit(1)
            print("Read bed file ...")
            read_feature = pd.read_table(bed, sep = "\t", header = 0)
            read_feature.columns = [c.strip("_") for c in read_feature.columns] 
            if len(set(read_feature["name"])) != len(read_feature):
                raise ValueError("Read name are not unqiue in the bed file.")
            read_feature["barcode"] = read_feature["name"].str.split(":", expand = True).iloc[:, -1]
            print("Collect cell level frame info ...")
            rna_df = pd.read_table(os.path.join(output_dir, sample + "_rna.stat"), sep = "\t", header = 0); rna_df.columns = ["barcode", "total_GFP"]
            # count reads for each barcode; 
            # count reads for each frame within each barcode
            # count read for each gene within each frame and barcode
            # pos for junction site
            read_feature["pos"] = np.where(read_feature["strand"] == "+", read_feature["start"], read_feature["end"])
            # count number of dominant fusion sites at single cell barcode for the same frame and same gene (one site required)
            barcode_pos_count = read_feature.groupby(['barcode', "frame", "gene"])["pos"].agg(mode_freq).reset_index(name = "pos_read")
            barcode_read = read_feature.groupby(["barcode"]).size().reset_index(name = "barcode_read"); frame_read = read_feature.groupby(["barcode", "frame"]).size().reset_index(name = "frame_read"); gene_read = read_feature.groupby(["barcode", "frame", "gene", "symbol"]).size().reset_index(name = "gene_read")
            # combine all information in one df
            barcode_stat = pd.merge(barcode_read, frame_read, on = "barcode")
            barcode_stat = pd.merge(barcode_stat, gene_read, on = ["barcode", "frame"])
            barcode_stat = pd.merge(barcode_stat, barcode_pos_count, on = ["barcode", "frame", "gene"])
            barcode_stat = pd.merge(rna_df, barcode_stat, on = "barcode", how = "outer")
            barcode_stat.to_csv(os.path.join(output_dir, sample + '_rna.gene.stat'), sep = "\t", header = True, index = False)
            # barcode_frame_stat["frame_ratio"] = round(barcode_frame_stat["frame_read"]/barcode_frame_stat["barcode_read"]*100)
            # barcode_frame_stat = pd.merge(barcode_stat, gene_read, on = ["barcode", "frame_"])
            # barcode_frame_stat["gene_frame_ratio"] = round(barcode_frame_stat["gene_frame_read"]/barcode_frame_stat["frame_read"]*100)
            # barcode_frame_stat.to_csv(os.path.join(output_dir, sample + ".gene.stat"), sep = "\t", header = True, index = False)
            # gene_clean = barcode_frame_stat.query("gene_frame_read >= 3 and gene_frame_ratio > 80 and frame_ratio > 80")
            # gene_clean.to_csv(os.path.join(output_dir, sample + ".gene.clean.stat"), sep = "\t", header = True, index = False)
            # target_frame_gene = gene_clean.query("frame_ == @frame")
            # print(f"Identified genes with target frame {frame}: {len(target_frame_gene)}")
            # print(f"target frame {frame} ratio: {round(len(target_frame_gene)/len(gene_clean)*100)}%")
            # print("Report genes ...")
            # from utilities.excel_tools import excel_write
            # excel_write(target_frame_gene, os.path.join(output_dir, sample + ".gene.xlsx"))
    if args.seq_type == "merge":
        print("Merge DNA & RNA stat ...")
        dna_df = pd.read_table(os.path.join(output_dir, sample + "_dna.stat"), sep = "\t", header = 0); dna_df.columns = ["ChiC_barcode", "ChiC_read"]
        rna_df = pd.read_table(os.path.join(output_dir, sample + "GFP_rna.gene.stat"), sep = '\t', header = 0)
        frame, frame_freq = np.unique(rna_df["frame"], return_counts=True)
        print(frame, frame_freq)
        dna_df["index"] = dna_df["ChiC_barcode"].str.split("+", expand = True)[0].astype(str) + "+" + dna_df["ChiC_barcode"].str.split("+", expand = True)[2].astype(str)
        rna_df["index"] = rna_df["barcode"].str.split("+", expand = True)[0].astype(str) + "+" + rna_df["barcode"].str.split("+", expand = True)[2].astype(str)
        dna_rna_combined = pd.merge(dna_df, rna_df, on = "index"); dna_rna_combined.drop("index", axis = 1, inplace = True)
        if args.cap is not None:
            print("Add CAP information ...")
            cap_df = pd.read_table(args.cap, sep = "\t", header = 0)
            dna_rna_combined["CAP"] = np.where(dna_rna_combined["symbol"].isin(cap_df["symbol"].tolist()), "CAP", "NULL")
        dna_rna_combined.to_csv(os.path.join(output_dir, sample + "_combined.stat"), sep = "\t", header = True, index = False)
        # filter 
        support_ratio = dna_rna_combined["pos_read"]/dna_rna_combined["barcode_read"]
        print(f"Apply frame as {args.frame} ...")
        dna_rna_filtered = dna_rna_combined[support_ratio > 0.7].query("frame == @args.frame and barcode_read > 30")
        dna_rna_filtered.to_csv(os.path.join(output_dir, sample + "_combined_filtered.stat"), sep = "\t", header = True, index = False)
        if args.cap is not None: 
            dna_rna_filtered.query("CAP == 'CAP'")[["ChiC_barcode", "symbol"]].to_csv(os.path.join(output_dir, sample + "_combined_cap.txt"), sep = "\t", header = False, index = False)

if __name__ == "__main__":
    main()
