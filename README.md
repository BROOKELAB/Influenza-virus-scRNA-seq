# Influenza-virus-scRNA-seq
Gene expression analysis of scRNA-seq libraries from A549s infected with H1N1 or H3N2.
### 10X Chromium single cell scripts

## Analysis Overview:
The scripts used in this analysis are specifically to look for differences in the percentage of
infected cells expressing selected genes of interest between H1N1 and H3N2. Scripts required for
reference assembly and overall gene expression analysis can be found in BROOKELAB/SingleCell. 

This analysis was performed using Cell Ranger on a 10X Chromium Single Cell experiment of influenza
infected human alveolar epithelial cells, A549. Some of the steps performed include quality
filtering, normalization, annotation, dimensional reduction,and quantification of gene expression
frequencies.

### Requirements:
The scripts requires R to be installed and made available from the command line * R (v4.1.2)

The following R packages are also required: * simpleSingleCell (v1.18.0) * DropletUtils (v1.14.2) *
scater (v1.22.0) * scran (v1.22.1) * sctransform (v0.3.5) * Seurat (v4.2.1) * ggplot2 (v3.4.0) *
dittoSeq (v1.6.0) * dplyr (v1.0.10) * glmGamPoi (v1.6.0)

### Preliminary steps: 
Raw reads should be demultiplexed and mapped to a host/virus hybrid reference using the 10X Chromium
Cell Ranger software package:
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

-1- Use cellranger counts for alignment, filtering, and barcode/UMI counting. It inputs FASTQ and
provides multiple file outputs (CSV, BAM, MEX, and H5). 
- txt file for input: files_*.txt
- Script for cellranger counts: cellranger_count_*.sh

-2- Perform cellranger aggr to aggregate all samples generated from one experiment, including
experimental conditions and replicates. This analysis requires .h5 files as input and provides a
list of tsv files and filtered matrixes that can be used for downstream analysis. 
- csv file for input: AggrList_all.csv 
- Script for cellranger aggr: aggr_all.sh


### Script: Cal07_Perth09_seurat.R

The script performs preliminary filtering (empty drops, cell cycle calling, filter cells by
detected genes, filter genes by min cells, double calling) and then calculates viral read
percentages to call infection status and viral gene presence/absence, resulting in normalized seurat
object for future analysis. It also generates expression plots of key antiviral genes (e.g. IFNL1
and IFIT3) to determine the level of expression in libraries.

The script uses the function calculate_expression_percentage to generat a CSV file containing
percentages of infected cells expressing diverse genes in Cal07 or Perth09 specific libraries. This
function needs to be run for each virus library. 
- Input for the function is a filtered seurat object containing only infected cells for either
		virus (cell_data) and a list of genes (gene_list). 
- Output is a csv file containing the percentage of cells expressing each gene from gene_list.
