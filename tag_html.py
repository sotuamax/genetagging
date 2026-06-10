#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os 
import pandas as pd
import argparse
import numpy as np 

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("-sample", "--sample", help="sample prefix used to find fastq file")
    parser.add_argument("-outdir", "--outdir", default = ".", help="output directory name")
    # parser.add_argument("-p", "--plot", help = "plot GFP read distribution", action = "store_true")
    # parser.add_argument("-cap", "--cap", help = "CAP file", required = False)
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    sample = args.sample
    output_dir = args.outdir
    # read GFP results 
    print("Read GFP results ...")
    rna_df = pd.read_table(os.path.join(output_dir, sample + "GFP_rna.gene.stat"), sep = "\t", header = 0)
    dna_df = pd.read_table(os.path.join(output_dir, sample + "_dna.stat"))
    combined_df = pd.read_table(os.path.join(output_dir, sample + "_combined.stat"), sep = "\t", header = 0)
    filtered_df = pd.read_table(os.path.join(output_dir, sample + "_combined_filtered.stat"), sep = "\t", header = 0)
    
    # about GFP_filter column: 
    # 1 -> ChiC w/ GFP support and GFP for single gene; 
    # 0 -> ChiC w/ GFP support and GFP not pass criteria; 
    # -1 -> ChiC wo/ GFP reads support.
    # freq_df = filtered_df.groupby(["symbol", "label"]).size().reset_index(name = "freq")
    # s_list = list()
    # for s, s_df in freq_df.groupby("symbol"):
    #     if len(s_df.query("label == 0")["freq"]) == 1:
    #         a = s_df.query("label == 0")["freq"].item() 
    #     else:
    #         a = 0
    #     if len(s_df.query("label == 1")["freq"]) == 1:
    #         b = s_df.query("label == 1")["freq"].item()
    #     else:
    #         b = 0
    #     s_list.append((s, a, b))
    # s_group = pd.DataFrame(s_list, columns = ["symbol", "groupA", "groupB"])
    # baseline = s_group["groupA"].sum()/s_group["groupB"].sum()
    # s_group["AB"] = (s_group["groupA"]/s_group["groupB"])/baseline
    # s_group["score"] = s_group["AB"]/baseline
    # gfp_df = pd.DataFrame({"Total reads":[round(total_GFP_read)], 
    #               "Valid primer (%)":[round(primer_yield,1)], 
    #               "Valid barcodes (%)":[round(barcode_yield, 1)], 
    #               "Valid reads (%)":[round(gfp_yiled, 1)]}).transpose()
    # gfp_gene_df = pd.DataFrame({"Valid frame (%)" : [target_frame], 
    #               "Total unique genes":[len(gene_freq)],
    #               "Top target genes":[", ".join([gene1, gene2, gene3])],
    #               "Top target freq":[", ".join([str(freq1), str(freq2), str(freq3)])]}).transpose()
    # read chic summary 
    total_ChiC_read = sum(chic_sublibrary_summary["total_reads"])
    ChiC_barcode_yield = chic_sublibrary_summary["barcode_yield"].mean()
    from utilities.parse_log import markdup_log_parser
    chic_dedup_df = pd.read_table(os.path.join(output_dir, sample + ".dedup.stat"), sep = "\t", header = 0, index_col = 0)
    duplication_ratio = chic_dedup_df.loc["PERCENT_DUPLICATION"].item()
    alignment_info = pd.read_table(os.path.join(output_dir, sample + ".flagstat"), header = None, names = ["pass", "fail", "category"])
    mapped_reads = alignment_info.query("category == 'primary mapped'")["pass"].item()
    mapped_perc = alignment_info.query("category == 'primary mapped %'")["pass"].item()
    # 
    dedup_info = pd.read_table(os.path.join(output_dir, sample + ".dedup.flagstat"), header = None, names = ["pass", "fail", "category"])
    dedup_mapped_reads = dedup_info.query("category == 'primary mapped'")["pass"].item()
    # 
    qc_info = pd.read_table(os.path.join(output_dir, sample + ".dedup.qc.flagstat"), header = None, names = ["pass", "fail", "category"])
    qc_dedup_mapped_reads = qc_info.query("category == 'primary mapped'")["pass"].item()
    chic_df = pd.DataFrame({"Total reads":[round(total_ChiC_read)], 
                            "Valid barcodes (%)":[round(ChiC_barcode_yield, 1)],
                            "Mapped reads" : [mapped_reads],
                            "Map ratio (%)":[mapped_perc],
                            "Reads after deduplication":[dedup_mapped_reads],
                            "Duplication (%)":[round(float(duplication_ratio)*100, 1)],
                            "Reads after deduplication+QC":[qc_dedup_mapped_reads]}).transpose()
    chic_barcode_distribution_png = os.path.join(output_dir, sample + "_barcode.png")
    import base64
    with open(chic_barcode_distribution_png, "rb") as img_file:
        encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    img_html = f'<img src="data:image/png;base64,{encoded_string}" width="800"/>'
    with open(os.path.join(output_dir, sample + "GFP.png"), "rb") as gfp_png:
        gfp_string = base64.b64encode(gfp_png.read()).decode('utf-8')
    gfp_img_html = f'<img src="data:image/png;base64,{gfp_string}" width="800"/>'
    # 
    with open(os.path.join(output_dir, sample + "_clean_barcode.png"), "rb") as img_file:
        encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    clean_img_html = f'<img src="data:image/png;base64,{encoded_string}" width="800"/>'
    # with open(os.path.join(output_dir, sample + "pair.png"), "rb") as img_file:
    #     encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    # pair_html = f'<img src="data:image/png;base64,{encoded_string}" width="800"/>'
    # read chic barcode stat 
    chic_barcode_stat = pd.read_table(os.path.join(output_dir, sample + "_barcode.stat"), sep = "\t", header = 0)
    chic_max = chic_barcode_stat["freq"].max(); chic_mean = chic_barcode_stat["freq"].mean(); chic_median = chic_barcode_stat["freq"].median()
    chic_img_stat = pd.DataFrame({"Total cell (no filtering)":[len(chic_barcode_stat)], 
                                  "Max reads per cell":[chic_max], 
                                  "Mean reads per cell":[chic_mean], 
                                  "Median reads per cell":[chic_median],
                                  }).transpose()
    chic_clean_img_stat = pd.DataFrame({"Total cells (filtered by GFP)":[len(GFP_ChiC_df)],
                                        "Valid cells (filtered by GFP & ChiC readcount)":[len(GFP_ChiC_df.query("label == 1"))],
                                        "Max reads per cell":[GFP_ChiC_df["ChiC_read_count"].max()], 
                                        "Mean reads per cell":[GFP_ChiC_df["ChiC_read_count"].mean()],
                                        "Median reads per cell":[GFP_ChiC_df["ChiC_read_count"].median()]}).transpose()
    #print(chic_img_stat)
    # tab pages
    gfp_summary = "<h3>GFP sequencing</h3>" + gfp_df.reset_index().to_html(index = False, header = False,border=0)
    gfp_gene_summary = "<h3>GFP gene</h3>" + gfp_gene_df.reset_index().to_html(index = False, header = False,border=0)
    gfp_plot = "<h3>GFP reads per cell</h3>" + gfp_img_html
    chic_summary = "<h3>ChiC sequencing</h3>" + chic_df.reset_index().to_html(index = False, header = False,border=0)
    img_summary = "<h3>ChiC reads per cell</h3>" + img_html
    clean_img_summary = "<h3>ChiC reads per cell filtered by GFP</h3>" + clean_img_html
    chic_img_stat_html = chic_img_stat.reset_index().to_html(index = False, header = False,border=0)
    chic_clean_img_stat_html = chic_clean_img_stat.reset_index().to_html(index = False, header = False, border = 0)
    # Build buttons/tabs (label)
    buttons_html = f"""
    <button class="tablinks" onclick="openTab(event, 'tab1')" id="defaultOpen1">Summary</button>
    <button class="tablinks" onclick="openTab(event, 'tab2')">GFP sublibrary summary</button>
    <button class="tablinks" onclick="openTab(event, 'tab3')">GFP tag genes</button>
    <button class="tablinks" onclick="openTab(event, 'tab4')">ChiC reads without any filtering</button>
    <button class="tablinks" onclick="openTab(event, 'tab5')">GFP tag genes frequency</button>
    <button class="tablinks" onclick="openTab(event, 'tab6')">GFP tag gene & ChiC read count</button>
    <button class="tablinks" onclick="openTab(event, 'tab7')">GFP tag gene group based on ChiC read</button>
    """
    gfp_sublibrary_summary_df.to_csv(os.path.join(output_dir, sample + "_gfp.csv"), sep = "\t", header = True, index = False)
    download_file = os.path.join(sample + "_gfp.csv")
    download_gfp = f"""
    <a href="{download_file}" download>
    <button style="margin: 10px 0;">Download CSV</button>
    </a>
    """
    gfp_tag_gene_df.to_csv(os.path.join(output_dir, sample + "_gene.csv") ,sep = "\t", header = True, index = False)
    download_file = os.path.join(sample + "_gene.csv")
    download_gene = f"""
    <a href="{download_file}" download>
    <button style="margin: 10px 0;">Download CSV</button>
    </a>
    """
    raw_ChiC_sum = raw_ChiC.groupby("GFP_filter")["ChiC_read"].sum().reset_index()
    raw_ChiC_sum["percentage(%)"] = round(raw_ChiC_sum["ChiC_read"]/sum(raw_ChiC_sum["ChiC_read"])*100, 2)
    raw_ChiC_sum["description"] = np.where(raw_ChiC_sum["GFP_filter"] == 1, "ChiC w/ GFP support for single gene", np.where(raw_ChiC_sum["GFP_filter"] == 0, "ChiC w/ GFP support but not for target frame or single gene", "ChiC wo/ GFP support"))
    raw_ChiC["GFP_filter"] = "GFP."+raw_ChiC["GFP_filter"].astype(str)
    ChiC_raw_tab = "<h3>ChiC reads w/wo GFP support</h3>" + raw_ChiC.to_html(index = False,border=0, table_id = "chic")
    ChiC_failed_GFP_tab = "<h3>ChiC failed GFP filter</h3>" + ChiC_failed_GFP.to_html(index = False,border=0, table_id = "chic_failed")
    # ChiC_raw_GFP0 = 
    ChiC_raw_summary  = "<h3>ChiC reads w/wo GFP summary</h3>" + raw_ChiC_sum.to_html(index = False, border=0)
    # Build contents html (label content)
    contents_html = f"""
    <div id="tab1" class="tabcontent">
    {gfp_summary}
    {gfp_gene_summary}
    {gfp_plot}
    {chic_summary}
    {img_summary}
    {chic_img_stat_html}
    {clean_img_summary}
    {chic_clean_img_stat_html}
    <br><br></div>
    <div id="tab2" class="tabcontent">{download_gfp}{gfp_sublibrary_summary_df.to_html(index = False,border=0, table_id = "gfp")}</div>
    <div id="tab3" class="tabcontent">{download_gene}{gfp_tag_gene_df.to_html(index = False,border=0, table_id = "gene")}</div>
    <div id="tab4" class="tabcontent">{ChiC_raw_summary}{ChiC_raw_tab}{ChiC_failed_GFP_tab}</div>
    <div id="tab5" class="tabcontent">{gene_freq.to_html(index = False,border=0, table_id = "freq")}</div>
    <div id="tab6" class="tabcontent">{GFP_ChiC_df.to_html(index = False, border = 0, table_id="combined")}</div>
    <div id="tab7" class="tabcontent">{s_group.to_html(index = False, border = 0, table_id="group")}</div>
    """
    # html template 
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="UTF-8" />
    <title>{sample}</title>
    <!-- DataTables CSS -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css" />

    <!-- jQuery -->
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>

    <!-- DataTables JS -->
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

    <style>
        .tab {{
        overflow: hidden;
        border-bottom: 1px solid #ccc;
        }}
        .tab button {{
        background-color: inherit;
        border: none;
        outline: none;
        padding: 10px 20px;
        cursor: pointer;
        transition: background-color 0.3s;
        font-size: 16px;
        }}
        .tab button:hover {{
        background-color: #ddd;
        }}
        .tab button.active {{
        background-color: #ccc;
        }}
        .tabcontent {{
        display: none;
        padding: 20px;
        border: 1px solid #ccc;
        border-top: none;
        }}
        table {{
        border-collapse: collapse;
        width: 50%;
        }}
        th, td {{
        border-bottom: 1px solid #ccc;  /* Row lines */
        padding: 8px;
        text-align: left;
        }}
        tr:last-child td {{
        border-bottom: none;  /* Optional: remove bottom border on last row */
        }}
    </style>
    </head>
    <body>

    <h2>Summary report for {sample}</h2>

    <div class="tab">
    {buttons_html}
    </div>

    {contents_html}
    
    <script>
    function openTab(evt, tabName) {{
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {{
        tabcontent[i].style.display = "none";
    }}
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {{
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }}
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
    if (tabName === 'tab2') {{
        if (!$.fn.DataTable.isDataTable('#gfp')) {{
        $('#gfp').DataTable();
        }}
    }}
    if (tabName === 'tab3') {{
        if (!$.fn.DataTable.isDataTable('#gene')) {{
        $('#gene').DataTable();
        }}
    }}
    if (tabName === 'tab4') {{
        if (!$.fn.DataTable.isDataTable('#chic')) {{
        $('#chic').DataTable();
        }}
        if (!$.fn.DataTable.isDataTable('#chic_failed')) {{
        $('#chic_failed').DataTable();
    }}
    }}
    if (tabName === 'tab5') {{
        if (!$.fn.DataTable.isDataTable('#freq')) {{
        $('#freq').DataTable();
        }}
    }}
    if (tabName === 'tab6') {{
        if (!$.fn.DataTable.isDataTable('#combined')) {{
        $('#combined').DataTable();
        }}
    }}
    if (tabName === 'tab7') {{
        if (!$.fn.DataTable.isDataTable('#group')) {{
        $('#group').DataTable();
        }}
    }}
    }}
    
    // Open the first tab by default
    document.getElementById("defaultOpen1").click();
    </script>

    </body>
    </html>
    """
    # write into html file 
    with open(os.path.join(output_dir, sample + ".html"), "w") as f:
        f.write(html_template)

if __name__ == "__main__":
    main()

