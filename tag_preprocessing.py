#!/usr/bin/env python3
import os 
import pandas as pd
import argparse
import logging
import glob 
import subprocess
"""
Preprocessing read by assigning primer and barcode:
primer carried by R1;
barcode carries by R2;
both primer/barcode are optional; 

Example:
mpiexec -n 12 tag_preprocessing.py -sample <GC> -dir <dir> -outdir <output> -primer [seq] -barcode [barcode]

For GFP, primer required; 
For ChiC, no primer given. 

"""

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("-sample", "--sample", help="sample prefix used to find fastq file")
    parser.add_argument("-outdir", "--outdir", default = ".", help="output directory name")
    parser.add_argument("-dir", "--directory", help = "read directory used to find its fastq file")
    parser.add_argument("-read_len", "--read_length", default = 100, type = int, help = "read length used for alignment")
    parser.add_argument("-primer", "--primer_seq", required = False, help = "primer sequence to be used for read filtering")
    # parser.add_argument("-regex", action = "store_true", help = "primer sequence provided in regular expression pattern.")
    parser.add_argument("-barcode", "--barcode", required = False, help = "cell barcode table used to filter and group on read")
    parser.add_argument("-p", "--prefix", action = "store_true", help = "if sample as prefix name (to add [0-9][0-9] for searching its R1/R2)")
    args=parser.parse_args()
    return args

def main():
    from mpi4py import MPI
    from mpi4py.util import pkl5
    comm = pkl5.Intracomm(MPI.COMM_WORLD) 
    # to overcome the overflow error when comm data > 2 GB
    rank = comm.Get_rank()
    size = comm.Get_size()
    args = args_parser()
    # requied argument
    output_dir = args.outdir
    try:
        os.mkdir(output_dir)
        print(f"Create {output_dir}")
    except Exception as e:
        print(e)
        pass
    sample = args.sample
    log_file = os.path.join(output_dir, sample +".log")
    rdir = args.directory
    if not args.prefix:
        r1_list = sorted(glob.glob(os.path.join(rdir, sample + "_R1.fastq.gz")))
        r2_list = sorted(glob.glob(os.path.join(rdir, sample + "_R2.fastq.gz")))
    else:
        r1_list = sorted(glob.glob(os.path.join(rdir, sample + "[0-9][0-9]_R1.fastq.gz")))
        r2_list = sorted(glob.glob(os.path.join(rdir, sample + "[0-9][0-9]_R2.fastq.gz")))
    if [r.split("_R1")[0] for r in r1_list] != [r.split("_R2")[0] for r in r2_list]:
        raise IOError("Cannot correctly locate paired R1/R2!")
    from utilities.fastq_tools import fq2df, read_len, read_select, write_read
    primer_seq = args.primer_seq.upper()
    if len(r1_list) == 0:
        raise IOError("No R1 data find")
    if rank == 0:
        print("Locate paired R1/R2 ...")
        w_idx = 0
        worker_tasks = {w:[] for w in range(size)}
        for (r1,r2) in zip(r1_list, r2_list):
            worker_tasks[w_idx].append((r1,r2))
            w_idx = (w_idx + 1) % size
        # with pre-processing step, start a new log file 
        logging.basicConfig(filename = log_file, filemode = "w", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.INFO)
        logging.info(f"Libraries\t{len(r1_list)}")
        print(f"Locate paired R1/R2 for the {sample} in {os.path.abspath(rdir)}")
        full_len = read_len(r1_list[0])
        logging.info(f"Raw read length\t{full_len}bp")
        pre_fw = open(os.path.join(output_dir, sample + ".pre.txt"), "a") 
        pre_fw.write("library\tindex\ttotal_reads\t")
        r1_head = 0; r1_tail = None
        if primer_seq is not None:
            pre_fw.write("primer_select\tprimer_yield\t")
            logging.info(f"Primer\t{primer_seq}")
            r1_head = len(primer_seq)
            logging.info(f"Primer length\t{r1_head}bp")
        barcode_list = list(); r2_head = 0
        if args.barcode is not None:
            barcode_list = pd.read_table(args.barcode, header = None)[0].tolist()
            r2_head = len(barcode_list[0])
            pre_fw.write("barcode_select\tbarcode_yield\t")
            logging.info(f"Barcode count\t{len(barcode_list)}")
        # write preprocessing log info 
        pre_fw.write("final_yield\n")
        if args.read_length + r1_head < full_len:
            r1_tail = r1_head + args.read_length
    else:
        worker_tasks = None
        r1_head = None
        r1_tail = None 
        r2_head = None 
        barcode_list = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    r1_head = comm.bcast(r1_head, root = 0)
    r1_tail = comm.bcast(r1_tail, root = 0)
    r2_head = comm.bcast(r2_head, root = 0)
    barcode_list = comm.bcast(barcode_list, root = 0)
    lib_list = list()
    for (r1,r2) in worker_tasks[rank]:
        lib_content = list()
        lib = os.path.basename(r1).split("_R1")[0]
        if os.path.exists(os.path.join(output_dir, lib + "_R1.fastq.gz")):
            print(f"Pass {lib} ...")
            continue
        print(f"Parse {lib} ... ")
        lib_content.append(lib)
        r1_df = fq2df(r1, quality = True)
        lib_label = ""
        lib_comment = sorted(set(r1_df["comment"]))
        if len(lib_comment) == 1:
            lib_label = lib_comment[0].split(":")[-1]
        lib_content.append(lib_label)
        total_read = len(r1_df)
        print(f"Total read {total_read}")
        lib_content.append(total_read)
        #pre_fw.write(f'{lib}\t{lib_label}\t{total_read}\t')
        if primer_seq is not None:
            print("Select on primer ...")
            r1_select_df = read_select(r1_df, [primer_seq])
            r1_select_read = len(r1_select_df)
            print(f"selected {r1_select_read} read")
            lib_content.append(r1_select_read)
            lib_content.append(round(r1_select_read/total_read*100))
            #pre_fw.write(f"{r1_select_read}\t{round(r1_select_read/total_read*100)}\t")
        else:
            r1_select_df = r1_df.copy()
        del r1_df 
        if args.barcode is not None:
            print("Select on barcode ...")
            r2_df = fq2df(r2, quality = True)
            r2_select_df = read_select(r2_df, barcode_list)
            r2_select_df["barcode"] = r2_select_df["seq"].str.slice(0, r2_head)
            r2_select_read = len(r2_select_df)
            lib_content.append(r2_select_read)
            lib_content.append(round(r2_select_read/total_read*100))
            # pre_fw.write(f"{r2_select_read}\t{round(r2_select_read/total_read*100)}\t")
            # get interaction
            select_read = list(set(r1_select_df["name"]).intersection(set(r2_select_df["name"])))
            r2_select_df = r2_select_df[r2_select_df["name"].isin(select_read)].copy()
            r2_select_df.sort_values(by = "name", inplace = True, ignore_index=True)
            r1_select_df = r1_select_df[r1_select_df["name"].isin(select_read)].copy()
            r1_select_df.sort_values(by = "name", inplace = True, ignore_index=True)
            if not r1_select_df['name'].equals(r2_select_df['name']):
                raise ValueError("R1 and R2 not match exactly on read name.")
            r1_select_df["barcode"] = r2_select_df["barcode"]
            r1_select_df["comment"] = r1_select_df["comment"].str.split(":", expand = True)[3]
            r1_select_df["comment"] = r1_select_df["comment"].astype(str) + "+" + r1_select_df["barcode"].astype(str)
            r2_select_df["comment"] = r2_select_df["comment"].str.split(":", expand = True)[3]
            r2_select_df["comment"] = r2_select_df["comment"].astype(str) + "+" + r2_select_df["barcode"].astype(str)
            r1_select_df.drop("barcode", axis = 1, inplace = True); r2_select_df.drop("barcode", axis = 1, inplace = True)
        else:
            r2_df = fq2df(r2, quality = True)
            r2_select_df = r2_df[r2_df["name"].isin(r1_select_df["name"])].copy()
        del r2_df
        # pre_fw.write(f"{round(len(r1_select_df)/total_read*100)}\n")
        lib_content.append(round(len(r1_select_df)/total_read*100))
        lib_list.append(lib_content)
        print("Write into new file ...")
        write_read(r1_select_df, output = os.path.join(output_dir, lib + "_R1.fastq.gz"), start=r1_head, end=r1_tail, comment = True)
        write_read(r2_select_df, output = os.path.join(output_dir, lib + "_R2.fastq.gz"), start=r2_head, comment = True)
    lib_list_all = comm.gather(lib_list, root = 0)
    ### 
    if rank == 0:
        for rank in lib_list_all:
            for lib in rank:
                cont = "\t".join([str(l) for l in lib])
                pre_fw.write(f"{cont}\n")
        if args.prefix:
            all_r1 = sorted(glob.glob(os.path.join(output_dir, sample + "[0-9][0-9]_R1.fastq.gz")))
            joint_r1 = " ".join(all_r1)
            all_r2 = sorted(glob.glob(os.path.join(output_dir, sample + "[0-9][0-9]_R2.fastq.gz")))
            joint_r2 = " ".join(all_r2)
            subprocess.call(f"cat {joint_r1} > {os.path.join(output_dir, sample)}_R1.fastq.gz", shell = True)
            subprocess.call(f"cat {joint_r2} > {os.path.join(output_dir, sample)}_R2.fastq.gz", shell = True)
            if os.path.exists(os.path.join(output_dir, sample + "_R1.fastq.gz")) and os.path.exists(os.path.join(output_dir, sample + "_R2.fastq.gz")):
                for rr1 in all_r1:
                    os.remove(rr1)
                for rr2 in all_r2:
                    os.remove(rr2)
            print(f"Combined R1/R2 from {len(r1_list)} libraries.")
            pre_fw.close()

if __name__ == "__main__":
    main()