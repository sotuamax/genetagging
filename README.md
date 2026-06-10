# gene tagging
leverage gRNA to tag target genes and fuse with GFP for the analysis of chromatin-associated protines. 

### Typical usage: 
- Preprocessing GFP fastq file to remove primer sequences, and label cell-barcode in title for following seq-alignment:
  - Use `tag_preprocessing.py -sample <GC>GFP -dir <file_directory> -outdir <output_directory> -primer CGTGGTCCTTATAGTCCATACC -barcode <known_cell_barcode.txt>`
- Perform seq-alignment for preprocessed GFP-seqs (see above):
  - Use: `tag_align.py -sample <GC>GFP -seq_type cDNA -outdir <output_directory> -ref <ref_genome.fa> -gtf <ref_genome.gtf> -n <threads>`
- Summarise GFP-target genes:
  - Use: `tag_summary.py -sample <GC>GFP -outdir <output_directory> -bc`
  - "-bc" enable cell-barcode mode; 
- Perform ChiC-seq reads alignment:
  - Use: `tag_align.py -sample <GC> -seq_type DNA -outdir <outdir> -ref <ref_genome.fa> -n <threads> -qc -report -GFP_filter`
  - "-qc" perform QC on alignment file;
  - "-report" generate summary report for ChiC-seq on single-cell basis;
- Report html file output:
  - Use `tag_html.py -sample <GC> -outdir <where_your_results_locate>`
    
