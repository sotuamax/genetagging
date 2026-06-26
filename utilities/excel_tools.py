import pandas as pd 

def read_sample(table):
    """
    Read sample table and filter for samples need to be run. 
    """
    if table.endswith(".txt"):
        sample_df = pd.read_table(table, sep = "\t", header = 0)
    if table.endswith(".xlsx"): 
        sample_df = pd.read_excel(table)
    return sample_df

def generate_alphabet(case:str):
    if case.lower() == "upper":
        return ''.join(chr(i) for i in range(ord('A'), ord('Z') + 1))
    if case.lower() == "lower":
        return ''.join(chr(i) for i in range(ord('a'), ord('z') + 1))

def excel_write(df, output, sheet_name = None):
    if not output.endswith(".xlsx"):
        output += ".xlsx"
    writer = pd.ExcelWriter(output, mode = "w")
    if sheet_name is not None:
        df.to_excel(writer, header = True, index = False, sheet_name = sheet_name)
    else:
        df.to_excel(writer, header = True, index = False)
    writer.close()

def excel_writer_style(sample_info, output):
    """
    Given all collected sample info, write into an excel file. 
    The newly created excel file contains: 
    sample_summary: all samples' read-statistics and gene number
    sample_gene: all samples' detected genes
    per sample gene: genes of all different frames detected per sample 
    """
    print("Write into excel ...")
    master_df, inframe_df, gene_df_dict  = sample_info
    # initiate excelWriter object 
    writer = pd.ExcelWriter(f"{output}.xlsx", mode = "w")
    # write master table into sample_summary sheet 
    master_df.to_excel(writer, sheet_name = "sample_summary", index = False)
    workbook = writer.book
    bold_format = workbook.add_format({'bold': True})
    worksheet = writer.sheets["sample_summary"]
    int_format = workbook.add_format({"num_format": "#,##0"})
    bg_format = workbook.add_format({"bg_color":"#fff7c7"})
    bg_bold_format = workbook.add_format({'bold': True, "bg_color":"orange"})
    worksheet.set_tab_color("red")
    for row in master_df.itertuples():
        if row.frame == "0":
            worksheet.write(row.Index+1, 6, master_df.iloc[row.Index, 6], bg_format)
            # worksheet.write(row.Index+1, 7, master_df.iloc[row.Index, 7], bg_format)
            if float(master_df.at[row.Index, "F0 read"].split("(")[-1].split("%")[0]) > 75:
                worksheet.write(row.Index+1, 7, master_df.iloc[row.Index, 7], bg_bold_format)
            else:
                worksheet.write(row.Index+1, 7, master_df.iloc[row.Index, 7], None)
        if row.frame == "1":
            worksheet.write(row.Index+1, 8, master_df.iloc[row.Index, 8], bg_format)
            if float(master_df.at[row.Index, "F1 read"].split("(")[-1].split("%")[0]) > 75:
                worksheet.write(row.Index+1, 9, master_df.iloc[row.Index, 9], bg_bold_format)
            else:
                worksheet.write(row.Index+1, 9, master_df.iloc[row.Index, 9], None)
        if row.frame == "2":
            worksheet.write(row.Index+1, 10, master_df.iloc[row.Index, 10], bg_format)
            if float(master_df.at[row.Index, "F2 read"].split("(")[-1].split("%")[0]) > 75:
                worksheet.write(row.Index+1, 11, master_df.iloc[row.Index, 11], bg_bold_format)
            else:
                worksheet.write(row.Index+1, 11, master_df.iloc[row.Index, 11], None)
    # set column format 
    worksheet.set_column(0, 1, 14, None)
    worksheet.set_column(1, 2, 16, int_format)
    worksheet.set_column(2, 12, 16, None)
    worksheet.set_column(12, 14, 10, None)
    # write sample genes into sample_gene sheet 
    inframe_df.to_excel(writer, sheet_name = "sample_gene", index = False)
    worksheet = writer.sheets["sample_gene"]
    # set column format 
    worksheet.set_column(0, 3, 14, None)
    pos_i = inframe_df.columns.tolist().index("pos")
    worksheet.set_column(pos_i, pos_i+1, 16, int_format)
    worksheet.set_column(pos_i+2, pos_i+3, 14, None)
    worksheet.set_tab_color("#FF9900")
    nrow,ncol = inframe_df.shape
    readcount_i = inframe_df.columns.tolist().index("read_count")
    readcount_col = generate_alphabet("upper")[readcount_i]
    worksheet.conditional_format(f"{readcount_col}1:{readcount_col}{nrow+1}", {'type':'cell', "criteria":">=", "value":5, "format":bold_format})
    unique_gene_total = len(set(inframe_df["symbol"]))
    unique_gene_cut5 = len(set(inframe_df.query("read_count >= 5")["symbol"]))
    unique_gene_cut10 = len(set(inframe_df.query("read_count >= 10")["symbol"]))
    unique_gene_cut20 = len(set(inframe_df.query("read_count >= 20")["symbol"]))
    unique_gene_cut50 = len(set(inframe_df.query("read_count >= 50")["symbol"]))
    unique_gene_cut100 = len(set(inframe_df.query("read_count >= 100")["symbol"]))
    anno_col = generate_alphabet("upper")[ncol + 2]
    worksheet.write(f"{anno_col}2", f"There are {unique_gene_total} unique genes;")
    worksheet.write(f"{anno_col}3", f"There are {unique_gene_cut5} unique genes with read-counts>=5;")
    worksheet.write(f"{anno_col}4", f"There are {unique_gene_cut10} unique genes with read-counts>=10;")
    worksheet.write(f"{anno_col}5", f"There are {unique_gene_cut20} unique genes with read-counts>=20;")
    worksheet.write(f"{anno_col}6", f"There are {unique_gene_cut50} unique genes with read-counts>=50;")
    worksheet.write(f"{anno_col}7", f"There are {unique_gene_cut100} unique genes with read-counts>=100;")
    # 
    unique_genes = inframe_df[["gene", "symbol", "freq"]].drop_duplicates()
    unique_genes.sort_values(by = ["freq", "symbol"], ascending = False, inplace = True, ignore_index = True)
    unique_genes.to_excel(writer, sheet_name = "unique_gene", index = False)
    worksheet = writer.sheets["unique_gene"]
    worksheet.set_column(1, 2, 16, None)
    worksheet.set_tab_color("yellow")
    # 
    for row in master_df.itertuples():
        sample_frame = int(row.frame)
        sample_gene = gene_df_dict[row.sample]
        sample_gene.to_excel(writer, sheet_name = row.sample, index = False)
        worksheet = writer.sheets[row.sample]
        pos_i = sample_gene.columns.tolist().index("pos")
        worksheet.set_column(pos_i, pos_i+1, 12, int_format)
        max_row,max_col = sample_gene.shape
        frame_i = sample_gene.columns.tolist().index("frame")
        frame_col = generate_alphabet("upper")[frame_i]
        worksheet.conditional_format(f"{frame_col}1:{frame_col}{max_row+1}", {'type':'cell', "criteria":"equal to", "value":int(sample_frame), "format":bg_format})
        sunique_gene_total = len(set(sample_gene["symbol"]))
        sunique_gene_cut5 = len(set(sample_gene.query("read_count >= 5 & frame == @sample_frame")["symbol"]))
        sunique_gene_cut10 = len(set(sample_gene.query("read_count >= 10 & frame == @sample_frame")["symbol"]))
        sunique_gene_cut20 = len(set(sample_gene.query("read_count >= 20 & frame == @sample_frame")["symbol"]))
        sunique_gene_cut50 = len(set(sample_gene.query("read_count >= 50 & frame == @sample_frame")["symbol"]))
        sunique_gene_cut100 = len(set(sample_gene.query("read_count >= 100 & frame == @sample_frame")["symbol"]))
        anno_col = generate_alphabet("upper")[max_col + 2]
        worksheet.write(f"{anno_col}2", f"There are {sunique_gene_total} unique genes inframe;")
        worksheet.write(f"{anno_col}3", f"There are {sunique_gene_cut5} unique genes inframe with read-counts>=5;")
        worksheet.write(f"{anno_col}4", f"There are {sunique_gene_cut10} unique genes inframe with read-counts>=10;")
        worksheet.write(f"{anno_col}5", f"There are {sunique_gene_cut20} unique genes inframe with read-counts>=20;")
        worksheet.write(f"{anno_col}6", f"There are {sunique_gene_cut50} unique genes inframe with read-counts>=50;")
        worksheet.write(f"{anno_col}7", f"There are {sunique_gene_cut100} unique genes inframe with read-counts>=100;")
    writer.close() 

