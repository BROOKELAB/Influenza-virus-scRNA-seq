#Codes to run on biocluster

# Bash codes:
# $ screen
# $ cd /home/labs/cbrooke_lab/Joel/sc_RNA_seq_2020/results
# $ srun -p normal --mem=100g -n 1 --pty bash
# $ module load  R/4.1.2-IGB-gcc-8.2.0
# $ R

#his starts interactive R session

library(Seurat)
library(sctransform)
library(glmGamPoi)
library(DropletUtils)
library(simpleSingleCell)
library(scater)
library(scran)
library(dplyr)


#Read aggregate file
all_data <- Read10X_h5("/outs/count/filtered_feature_bc_matrix.h5")

# Make Seurat object. The names.field = 2, names.delim = "-" will split off the number from the cellbarcode and use it for orig.ident
allcells_prefilt <- CreateSeuratObject(counts = all_data, min.cells = 4,
                                       min.features = 400,
                                       project = "Cal07_Perth09", names.field = 2, names.delim = "-")

# Read in the aggregation csv
sample_order <- read.csv("/outs/aggregation.csv")

# Change sample info
sample_order$sample_id[1:16] <- c("Cal07_bystander_16hr_1","Cal07_bystander_16hr_2","Cal07_bystander_8hr_1","Cal07_bystander_8hr_2",
                                  "Cal07_infected_16hr_1","Cal07_infected_16hr_2","Cal07_infected_8hr","Cal07_mock","Perth09_bystander_16hr_1","Perth09_bystander_16hr_2",
                                  "Perth09_bystander_8hr_1","Perth09_bystander_8hr_2","Perth09_infected_16hr_1","Perth09_infected_16hr_2","Perth09_infected_8hr",
                                  "Perth09_mock")

# The “numbers” in the Seurat object are treated as character so once you have >9 samples, they won’t be in the correct order. Instead do the same in the sampOrder:
sample_order$number <- as.character(1:nrow(sample_order))

# orig.ident is a factor, so it’s easy to change the levels of the factor to the sample_id:
levels(allcells_prefilt$orig.ident) <- sample_order$sample_id[order(sample_order$number)]

# change the order of the levels to the order in sample_order (assuming that the order you want)
allcells_prefilt$orig.ident <- factor(allcells_prefilt$orig.ident, levels = sample_order$sample_id )

# add in the other info columns that you had in your AggrList_all.csv:

rownames(sample_order) <- sample_order$sample_id
allcells_prefilt $group <- sample_order[allcells_prefilt $orig.ident, "group"]
allcells_prefilt $time <- sample_order[allcells_prefilt $orig.ident, "time"]
allcells_prefilt $rep <- sample_order[allcells_prefilt $orig.ident, "rep"]
allcells_prefilt $virus <- sample_order[allcells_prefilt $orig.ident, "virus"]

#QC

allcells_prefilt[["percent.mt"]] <- PercentageFeatureSet(allcells_prefilt, pattern = "^MT-")

##Remove mitochondria genes to get rid of dead cells

allcells_prefilt <- subset(allcells_prefilt, subset = percent.mt < 25)

VlnPlot(allcells_prefilt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("data_mtQC_101623.pdf")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(allcells_prefilt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(allcells_prefilt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("data_featureQC_101623.pdf",plot = plot1 + plot2)

#Save before cell cycle

saveRDS(allcells_prefilt,"101623_data_seurat_preCellCycle.rds")

#Read files

allcells_prefilt <- readRDS("101623_data_seurat_preCellCycle.rds")

#cell cycle genes

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

allcells_prefilt <- CellCycleScoring(allcells_prefilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#view cell cycle scores

head(allcells_prefilt[[]])

#Save copy of data post all QC

saveRDS(allcells_prefilt, "101623_allcells_postQC.rds")

#Read file to start Virus count calling

allcells_postfilt <- readRDS("101623_allcells_postQC.rds")

# viral gene names

segment_gene_names <- c("NS", "NA.", "HA", "NP", "PB2", "PB1", "PA", "M")

#### Calculate host counts and viral counts

temp <- GetAssayData(allcells_postfilt, slot = "count")

allcells_postfilt$host_counts <- colSums(temp[!rownames(temp) %in% segment_gene_names,])
allcells_postfilt$virus_counts <- colSums(temp[segment_gene_names,])
allcells_postfilt$virus_pct <- allcells_postfilt$virus_counts / allcells_postfilt$nCount_RNA *100

#look at distributions on log10 scale
#get half the minimum non-zero value

halfmin_virus <- min(allcells_postfilt$virus_pct[allcells_postfilt$virus_pct>0])/2
allcells_postfilt$virus_pct_log10 <- log10(allcells_postfilt$virus_pct+halfmin_virus)

#Try to find local minima points
des.all <- density(allcells_postfilt$virus_pct_log10[allcells_postfilt$virus_counts > 0])
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]

Idents(allcells_postfilt) <- "All"
p <- RidgePlot(allcells_postfilt[,allcells_postfilt$virus_counts > 1], features = "virus_pct_log10",  assay = "RNA",
               same.y.lims = TRUE,
               log = FALSE, slot = "data")

ggsave("virus_pct_log10_RidgePlot_Perth_timecourse_withMinima_101623.jpeg",
       plot = p, width = 8, height = 5, units = "in", type = "cairo")

#there are multiple dips; based on visual inspection, use
min.all.use <- -1.25

#call infected any cell above the first minimum ----

allcells_postfilt$InfectedStatus <- "NotInfected"
allcells_postfilt$InfectedStatus[allcells_postfilt$virus_pct_log10 > min.all.use] <- "Infected"

#Move viral genes to cell metadata

temp <- GetAssayData(allcells_postfilt, slot = "count")[segment_gene_names,] %>% as.matrix() %>% t()  %>% as.data.frame()

allcells_postfilt <- AddMetaData(allcells_postfilt, temp)

allcells_postfilt <- allcells_postfilt[!rownames(allcells_postfilt) %in% segment_gene_names]

#Normalization

allcells_postfilt <- NormalizeData(allcells_postfilt, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable gene

allcells_postfilt <- FindVariableFeatures(allcells_postfilt, selection.method = "vst", nfeatures = 2000)


#Scale Data after normalization

allcells_postfilt <- SCTransform(allcells_postfilt,
                                 method = "glmGamPoi",
                                 vars.to.regress = c("S.Score", "G2M.Score"),
                                 return.only.var.genes = TRUE)

#Run pca

allcells_postfilt <- RunPCA(allcells_postfilt, verbose = FALSE)

print(allcells_postfilt[["pca"]], dims = 1:5, nfeatures = 5)

#findn cluster

allcells_postfilt <- FindNeighbors(allcells_postfilt, dims = 1:40)
allcells_postfilt <- FindClusters(allcells_postfilt, resolution = 0.3)


#Do UMAP showing 0.3 res clusters

allcells_postfilt <- RunUMAP(allcells_postfilt,
                             dims = 1:40,
                             verbose = FALSE)

DimPlot(allcells_postfilt, reduction = "umap")
ggsave("data_UMAP_res0.3_101823.pdf")


#Save post UMAP

saveRDS(allcells_postfilt, "data_seurat_postUMAP_101823.rds")

#read data

allcells_postfilt <- readRDS("data_seurat_postUMAP_101823.rds")


#split by virus

DimPlot(allcells_postfilt, split.by = 'virus', group.by = 'InfectedStatus')

ggsave("split_UMAP_virus_101823.pdf")


#look at expression of IFNL1

IFNL1_cells <- WhichCells(allcells_postfilt, expression = IFNL1>0)

DimPlot(allcells_postfilt, split.by = 'virus', cells.highlight= IFNL1_cells, cols.highlight = 'red', cols= "grey", pt.size = 2)

ggsave("split_virus_IFNL1_101124.png", width = 30, height = 10, units = "cm", dpi = 600, device = "png", type = "cairo")

#NP/NS and IFNL1 expression

expression_counts <- FetchData(allcells_postfilt, vars = c("IFNL1","NP"), slot = "counts")

expression_counts$virus <- allcells_postfilt@meta.data$virus

NP_INFL1_virus <- ggplot(expression_counts, aes(x = IFNL1, y = NP, color = virus)) +
  geom_point(size = 1) +
  xlab("log10(IFNL1 Counts)") +
  ylab("log10(NP Counts)") +
  theme_classic() +
  scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  scale_x_log10() +  # Set x-axis limits and breaks
  scale_y_log10()    # Set y-axis limits and breaks

ggsave("IFNL1_NP_101124.svg",plot = NP_INFL1_virus, width = 8, height = 5, units = "cm", dpi = 600, device = "svg")

#look at expression of IFIT3

FeaturePlot(allcells_postfilt, features = 'IFIT3', split.by ='group')

ggsave("IFIT3_split.pdf",width = 30, height = 10, units = "cm")

IFIT3_cells <- subset(allcells_postfilt, subset = IFIT3 >0, slot = "counts")

table(IFIT3_cells$group)

IFIT3_cells <- WhichCells(allcells_postfilt, expression = IFIT3>0)

DimPlot(allcells_postfilt, split.by = 'virus', cells.highlight= IFIT3_cells, cols.highlight = 'red', cols= "grey")

ggsave("split_virus_IFIT3_052124.pdf",width = 30, height = 10, units = "cm")


#calculate % of infected cells expressing genes of interest

# Create function

calculate_expression_percentage <- function(cell_data, gene_list) {
  
  # Get the variable names 
  available_genes <- rownames(cell_data)
  
  # Create empty data frame to store the results
  results_df <- data.frame(Gene = character(), Percentage = numeric(), stringsAsFactors = FALSE)
  
  # Integrate the genes from gene_list
  for (gene in gene_list) {
    
    # Find the gene (case insensitive)
    matched_gene <- available_genes[grep(paste0("^", gene, "$"), available_genes, ignore.case = TRUE)]
    
    # Check if the gene was found
    if (length(matched_gene) == 1) {
      
      # Fetch data for the given gene
      gene_data <- FetchData(cell_data, vars = matched_gene, slot = "counts")
      
      # Calculate the percentage of cells expressing the gene in the infected group
      infected_cells <- cell_data$group == "infected"
      expressing_cells <- gene_data[infected_cells, matched_gene] > 0
      expression_percentage <- sum(expressing_cells) / sum(infected_cells) * 100
      
      # Add the result to the data frame
      results_df <- rbind(results_df, data.frame(Gene = gene, Percentage = expression_percentage))
    } else {
      
      # Append the result as NA if the gene is not found or not uniquely matched
      results_df <- rbind(results_df, data.frame(Gene = gene, Percentage = NA))
    }
  }
  
  # Return the results data frame
  return(results_df)
}

# Subset the cells for the specific virus
Cal07_cells <- subset(allcells_postfilt, subset = virus == "Cal07" & time == "16hr")
Perth09_cells <- subset(allcells_postfilt, subset = virus == "Perth09" & time == "16hr")

# Define the list of genes
gene_list <- read.csv("list_isgs.csv")
gene_list <- as.vector(gene_list$ISGs)

# Call the function and save
expression_percentages_df <- calculate_expression_percentage(Cal07_cells, gene_list)
write.csv(expression_percentages_df, file = "cal07_expression_percentages.csv", row.names = FALSE)

expression_percentages_df <- calculate_expression_percentage(Perth09_cells, gene_list)
write.csv(expression_percentages_df, file = "perth09_expression_percentages.csv", row.names = FALSE)

#Correlation NS1 and IFNL1/IFIT3

data_corr <-read.csv("H3N2_NS1_IFN.csv")

cor_test <- cor.test(data_corr$IFNL1, data_corr$NS1, method = "pearson")
correlation <- cor_test$estimate
p_value <- cor_test$p.value

corr_plot <- ggpubr::ggscatter(data_corr, x = "NS1", y = "IFNL1", 
                          add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                          conf.int = TRUE, cor.coef = TRUE, cor.method = 'pearson', size = 1,
                          xlab = 'Ratio NS1/NP Normalized to ACTB', ylab = 'IFNL1 Fold Change Relative to 2009', 
                          cor.coeff.args = list(method = "pearson", label.x.npc = "middle", label.y.npc = "top"))

ggsave("correlation_H3_IFNL1_101124.svg",plot = corr_plot, width = 8, height = 5, units = "cm", dpi = 600, device = "svg",
       path = "")

sessionInfo()

R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Linux (unknown distro)

Matrix products: default
BLAS/LAPACK: /home/apps/software/OpenBLAS/0.3.5-GCC-8.2.0-2.32/lib/libopenblas_sandybridgep-r0.3.5.so

locale:
  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
[9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
  [1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
  [1] cowplot_1.1.1               edgeR_3.36.0               
[3] limma_3.50.0                ComplexHeatmap_2.10.0      
[5] magrittr_2.0.3              BiocSingular_1.10.0        
[7] scran_1.22.1                scater_1.22.0              
[9] scuttle_1.4.0               simpleSingleCell_1.18.0    
[11] DropletUtils_1.14.2         SingleCellExperiment_1.16.0
[13] SummarizedExperiment_1.24.0 Biobase_2.54.0             
[15] GenomicRanges_1.46.1        GenomeInfoDb_1.30.0        
[17] IRanges_2.28.0              S4Vectors_0.32.3           
[19] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
[21] matrixStats_0.62.0          glmGamPoi_1.6.0            
[23] sctransform_0.3.5           SeuratObject_4.1.3         
[25] Seurat_4.2.1                forcats_0.5.2              
[27] stringr_1.4.1               dplyr_1.0.10               
[29] purrr_0.3.5                 readr_2.1.3                
[31] tidyr_1.2.1                 tibble_3.1.8               
[33] tidyverse_1.3.2             dittoSeq_1.6.0             
[35] ggplot2_3.4.0              

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                spatstat.explore_3.0-3   
[3] reticulate_1.26           R.utils_2.12.1           
[5] tidyselect_1.2.0          htmlwidgets_1.5.4        
[7] BiocParallel_1.28.2       Rtsne_0.16               
[9] ScaledMatrix_1.2.0        munsell_0.5.0            
[11] codetools_0.2-18          ica_1.0-3                
[13] statmod_1.4.37            future_1.29.0            
[15] miniUI_0.1.1.1            withr_2.5.0              
[17] spatstat.random_3.0-1     colorspace_2.0-3         
[19] progressr_0.11.0          knitr_1.40               
[21] ROCR_1.0-11               tensor_1.5               
[23] listenv_0.8.0             GenomeInfoDbData_1.2.7   
[25] polyclip_1.10-4           pheatmap_1.0.12          
[27] rhdf5_2.38.0              parallelly_1.32.1        
[29] vctrs_0.5.0               generics_0.1.3           
[31] xfun_0.41                 timechange_0.1.1         
[33] doParallel_1.0.17         R6_2.5.1                 
[35] clue_0.3-65               ggbeeswarm_0.7.2         
[37] rsvd_1.0.5                locfit_1.5-9.6           
[39] bitops_1.0-7              rhdf5filters_1.6.0       
[41] spatstat.utils_3.0-1      DelayedArray_0.20.0      
[43] assertthat_0.2.1          promises_1.2.0.1         
[45] scales_1.2.1              googlesheets4_1.0.1      
[47] beeswarm_0.4.0            gtable_0.3.1             
[49] beachmat_2.10.0           globals_0.16.1           
[51] processx_3.8.0            goftest_1.2-3            
[53] rlang_1.0.6               GlobalOptions_0.1.2      
[55] splines_4.1.2             lazyeval_0.2.2           
[57] gargle_1.2.1              spatstat.geom_3.0-3      
[59] broom_1.0.1               reshape2_1.4.4           
[61] abind_1.4-5               modelr_0.1.9             
[63] backports_1.4.1           httpuv_1.6.6             
[65] tools_4.1.2               ellipsis_0.3.2           
[67] RColorBrewer_1.1-3        ggridges_0.5.4           
[69] Rcpp_1.0.9                plyr_1.8.7               
[71] sparseMatrixStats_1.6.0   zlibbioc_1.40.0          
[73] RCurl_1.98-1.9            ps_1.7.2                 
[75] deldir_1.0-6              GetoptLong_1.0.5         
[77] viridis_0.6.2             pbapply_1.5-0            
[79] zoo_1.8-11                haven_2.5.1              
[81] ggrepel_0.9.2             cluster_2.1.4            
[83] fs_1.5.2                  data.table_1.14.4        
[85] scattermore_0.8           circlize_0.4.15          
[87] lmtest_0.9-40             reprex_2.0.2             
[89] RANN_2.6.1                googledrive_2.0.0        
[91] fitdistrplus_1.1-8        evaluate_0.18            
[93] hms_1.1.2                 patchwork_1.1.2          
[95] mime_0.12                 xtable_1.8-4             
[97] XML_3.99-0.12             readxl_1.4.1             
[99] shape_1.4.6               gridExtra_2.3            
[101] compiler_4.1.2            KernSmooth_2.23-20       
[103] crayon_1.5.2              R.oo_1.25.0              
[105] htmltools_0.5.3           later_1.3.0              
[107] tzdb_0.3.0                lubridate_1.9.0          
[109] DBI_1.1.3                 dbplyr_2.2.1             
[111] MASS_7.3-58.1             Matrix_1.5-1             
[113] cli_3.4.1                 R.methodsS3_1.8.2        
[115] metapod_1.2.0             parallel_4.1.2           
[117] igraph_1.3.5              pkgconfig_2.0.3          
[119] sp_1.5-1                  plotly_4.10.1            
[121] spatstat.sparse_3.0-0     foreach_1.5.2            
[123] xml2_1.3.3                vipor_0.4.5              
[125] dqrng_0.3.0               XVector_0.34.0           
[127] rvest_1.0.3               callr_3.7.3              
[129] digest_0.6.30             RcppAnnoy_0.0.20         
[131] graph_1.72.0              spatstat.data_3.0-0      
[133] rmarkdown_2.18            cellranger_1.1.0         
[135] leiden_0.4.3              uwot_0.1.14              
[137] DelayedMatrixStats_1.16.0 shiny_1.7.3              
[139] rjson_0.2.21              lifecycle_1.0.3          
[141] nlme_3.1-160              jsonlite_1.8.3           
[143] Rhdf5lib_1.16.0           BiocNeighbors_1.12.0     
[145] CodeDepends_0.6.5         viridisLite_0.4.1        
[147] fansi_1.0.3               pillar_1.8.1             
[149] lattice_0.20-45           fastmap_1.1.0            
[151] httr_1.4.4                survival_3.4-0           
[153] glue_1.6.2                iterators_1.0.14         
[155] png_0.1-7                 bluster_1.4.0            
[157] stringi_1.7.8             HDF5Array_1.22.1         
[159] irlba_2.3.5.1             future.apply_1.10.0      


