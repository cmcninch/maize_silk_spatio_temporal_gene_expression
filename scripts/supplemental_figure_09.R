### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(readxl)
library(RColorBrewer)
library(pheatmap)

### Load data ###
#################################################################################################################
library_codes <- fread("./data/Library_Codes.csv") %>%
  mutate_at(vars(-Library_Name), factor)

b73_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv")  %>%
  dplyr::select(GeneID, (filter(library_codes, Genotype == "B73"))$Library_Name)

mo17_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, (filter(library_codes, Genotype == "Mo17"))$Library_Name)

all_rpkm <- bind_rows(b73_rpkm, mo17_rpkm) %>% gather("Library_Name", "RPKM", 2:81)
all_rpkm <- left_join(all_rpkm, library_codes, by = "Library_Name")

hormone_genes <- bind_rows((fread("./data/Hormone_Biosynthesis_and_Marker_Genes.csv") %>%
                                dplyr::select(`Gene Name`, Order, `B73 GeneID`) %>% dplyr::rename("GeneID" = `B73 GeneID`) %>% mutate(Genotype = "B73")), 
                             (fread("./data/Hormone_Biosynthesis_and_Marker_Genes.csv") %>%
                                dplyr::select(`Gene Name`, Order, `Mo17 Best Blast Hit`) %>% dplyr::rename("GeneID" = `Mo17 Best Blast Hit`) %>% mutate(Genotype = "Mo17")))

### Define functions ###
#################################################################################################################
# Differential expression analysis
DiffExp_func <- function(count_data, genotypes = c("B73", "Mo17"), sections = c("A", "B", "C", "D", "E"), silk_locations = c("Emerged", "Husk-Encased"), years = c("2014", "2015"), 
                         design, contrast, read_depth_cutoff = 5, prop_libs_with_reads = 0.25, FDR = 0.05, FC = 2){
  
  sample_info <- filter(library_codes, UQ(as.name(contrast[1])) %in% contrast[2:3], Genotype %in% genotypes, Section %in% sections, Year %in% years, Location_Of_Silk %in% silk_locations) %>%
    arrange(UQ(as.name(contrast[1]))) %>% column_to_rownames(var = "Library_Name")
  count_data <- dplyr::select(count_data, GeneID, rownames(sample_info)) %>% arrange(GeneID)
  max_low_expressed_libraries <- (1 - prop_libs_with_reads) * (ncol(count_data) -1)
  disqualified_genes <- count_data[rowSums(count_data < read_depth_cutoff) > max_low_expressed_libraries, ]$GeneID
  count_data <- count_data %>% column_to_rownames(var = "GeneID")
  ddsDF <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = design) 
  ddsDF <- DESeq(ddsDF)
  `%notIN%` <- Negate(`%in%`)
  res <- as.data.frame(results(ddsDF, contrast = contrast, alpha = FDR)) %>% rownames_to_column(var = "GeneID") %>% filter(padj <= FDR, GeneID %notIN% disqualified_genes)
  res <- bind_rows(filter(res, log2FoldChange >= log2(FC)), filter(res, log2FoldChange <= -log2(FC))) %>%
    dplyr::select(GeneID:log2FoldChange, padj) %>%
    mutate(Expression_Change = ifelse(log2FoldChange > 0, "Expression Decreases", "Expression Increases"))
  return(res)
}

### Compute stats ###
#################################################################################################################
hormone_genes_exp <- bind_rows(mutate_at(b73_rpkm, vars(-GeneID), .funs = funs(log2((. + 1)))) %>%
                                   filter(GeneID %in% filter(hormone_genes, Genotype == "B73")$GeneID) %>%
                                   gather("Library_Name", "Log2RPKM", 2:41) %>%
                                   inner_join(., dplyr::select(library_codes, Library_Name, Section), by = "Library_Name") %>%
                                   group_by(GeneID, Section) %>% summarise(Log2RPKM = mean(Log2RPKM, na.rm = TRUE)),
                                 mutate_at(mo17_rpkm, vars(-GeneID), .funs = funs(log2((. + 1)))) %>%
                                   filter(GeneID %in% filter(hormone_genes, Genotype == "Mo17")$GeneID) %>%
                                   gather("Library_Name", "Log2RPKM", 2:41) %>%
                                   inner_join(., dplyr::select(library_codes, Library_Name, Section), by = "Library_Name") %>%
                                   group_by(GeneID, Section) %>% summarise(Log2RPKM = mean(Log2RPKM, na.rm = TRUE))) %>%
  ungroup() %>% inner_join(., hormone_genes, by = "GeneID") %>% 
  mutate(Section = paste(Genotype, "_", Section, sep = "")) %>% dplyr::select(`Gene Name`, Order, Section, Log2RPKM) %>%
  spread(Section, Log2RPKM) %>% ungroup() %>% arrange(Order)

hormone_genes_exp <- as.data.frame(hormone_genes_exp %>% filter(Order != "NA") %>% 
                                      dplyr::select(-Order) %>% column_to_rownames(var = "Gene Name"))

### Graph Stats ###
#################################################################################################################
colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)

colfunc <- colorRampPalette(c("yellow", "firebrick"))

pheatmap(hormone_genes_exp[, 1:5],
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         gaps_row = c(8, 13, 27),
         gaps_col = 5,
         scale = "row",
         cellwidth = 30,
         cellheight = 20,
         show_rownames = FALSE, 
         show_colnames = FALSE,
         col = colfunc(15),
         legend = FALSE,
         legend_labels = NULL,
         border_color = "black")

pheatmap(hormone_genes_exp[, 6:10],
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         gaps_row = c(8, 13, 27),
         gaps_col = 5,
         scale = "row",
         cellwidth = 30,
         cellheight = 20,
         show_rownames = FALSE, 
         show_colnames = FALSE,
         col = colfunc(15),
         legend = FALSE,
         legend_labels = NULL,
         border_color = "black")