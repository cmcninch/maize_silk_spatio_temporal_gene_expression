### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(DESeq2)
library(stats)

### Load data ###
#################################################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
  mutate_at(vars(-Library_Name), factor)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name)

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% 
  distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% 
  mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

### Define functions ###
#################################################################################################################
DiffExp_func <- function(count_data, genotypes = c("B73", "Mo17"), sections = c("A", "B", "C", "D", "E"), silk_locations = c("Emerged", "Husk-Encased"), years = c("2014", "2015"), 
                         harvest_date = c("7/30/14", "7/30/15", "8/2/15", "7/31/15", "8/4/15"), design, contrast, read_depth_cutoff = 5, prop_libs_with_reads = 0.25){
  
  sample_info <- filter(library_info, UQ(as.name(contrast[1])) %in% contrast[2:3], Genotype %in% genotypes, Section %in% sections, Year %in% years, Location_Of_Silk %in% silk_locations, 
                        Harvest_Date %in% Harvest_Date) %>%
    arrange(UQ(as.name(contrast[1]))) %>% column_to_rownames(var = "Library_Name")
  count_data <- dplyr::select(count_data, GeneID, rownames(sample_info)) %>% arrange(GeneID)
  max_low_expressed_libraries <- (1 - prop_libs_with_reads) * (ncol(count_data) -1)
  count_data <- count_data %>% column_to_rownames(var = "GeneID")
  ddsDF <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = design) 
  ddsDF <- DESeq(ddsDF)
  `%notIN%` <- Negate(`%in%`)
  res <- as.data.frame(results(ddsDF, contrast = contrast, alpha = 0.999)) %>% rownames_to_column(var = "GeneID")
  return(res)
}

scatter_plot_func <- function(df, box_df, draw_boxes, cor_value, p_value){
  
  plot <- ggplot(df, aes(x = B73_Log2FC, y = Mo17_Log2FC)) +
          scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1), expand = c(0, 0)) +
          scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1), expand = c(0, 0)) +
          theme_bw() +
          theme(axis.text = element_text(size = 24), 
                axis.title.x = element_text(size = 32,  margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
                panel.border = element_rect(size = 3), 
                legend.text = element_text(size = 16), 
                legend.title = element_text(size = 20),
                plot.margin = unit(c(1.2, 1, 1, 1),"cm")) +
          labs(x = "B73 Log2 Fold Change", 
               y = "Mo17 Log2 Fold Change", 
               fill = NULL)
  
  if (draw_boxes == FALSE) {
    print(plot + geom_point(alpha = 0.5)  +
            annotate(x = -3, y = 3.75, geom = "text", 
                     label = paste("r = ", cor_value, ";\np < ", p_value, sep = ""), size = 12))
    } else {
    print(plot + geom_rect(data = box_df, 
                           mapping = aes(xmin = xlow, xmax = xhigh, ymin = ylow, ymax = yhigh, fill = Type), 
                           alpha = 0.4, color = "black") +
            geom_point(alpha = 0.5) +
            annotate(x = -3, y = 3.75, geom = "text", 
                     label = paste("r = ", cor_value, ";\np < ", p_value, sep = ""), size = 12))
    }
}

### Compute Stats ###
#################################################################################################################
# Differential expression analysis at c to d silk transition
b73_c_to_d <- DiffExp_func(count_data = b73_counts, genotypes = "B73", 
                           design = ~ Section + Year, contrast = c("Section", "C", "D")) %>%
  dplyr::select(GeneID, log2FoldChange, padj) %>%
  dplyr::rename(B73_V4_GeneID = GeneID, 
                B73_Log2FC = log2FoldChange, 
                B73_FDR = padj) %>%
  left_join(., syntenic_gene_pairs, by = "B73_V4_GeneID")

mo17_c_to_d <- DiffExp_func(count_data = mo17_counts, genotypes = "Mo17", 
                            design = ~ Section + Year, contrast = c("Section", "C", "D")) %>%
  dplyr::select(GeneID, log2FoldChange, padj) %>%
  dplyr::rename(Mo17_CAU_GeneID = GeneID, 
                Mo17_Log2FC = log2FoldChange, 
                Mo17_FDR = padj) %>%
  left_join(., syntenic_gene_pairs, by = "Mo17_CAU_GeneID")

c_to_d <- left_join(dplyr::select(b73_c_to_d, Syntenic_Pair, B73_Log2FC, B73_FDR) %>% drop_na(), 
                    dplyr::select(mo17_c_to_d, Syntenic_Pair, Mo17_Log2FC, Mo17_FDR) %>% drop_na(), by = "Syntenic_Pair") %>%
  drop_na() %>%
  mutate(B73_Log2FC = -1*B73_Log2FC, 
         Mo17_Log2FC = -1*Mo17_Log2FC)

boxes <- data.frame(xlow = c(-1, -1, 1, -5, -5, 1, -5, 1, -1), 
                    xhigh = c(1, 1, 5, -1, -1, 5, -1, 5, 1), 
                    ylow = c(1, -5, -1, -1, -5, 1, 1, -5, -1), 
                    yhigh = c(5, -1, 1, 1, -1, 5, 5, -1, 1),
                    Type = c("Mo17 specific differential expression", "Mo17 specific differential expression",
                             "B73 specific differential expression", "B73 specific differential expression", 
                             "Concordant differential expression", "Concordant differential expression", 
                             "Discordant differential expression", "Discordant differential expression", 
                             "No differential expression"),
                    B73_Log2FC = 1,
                    Mo17_Log2FC = 1) %>%
  mutate(Type = factor(Type, 
                       levels = c("B73 specific differential expression", "Mo17 specific differential expression",
                                  "Concordant differential expression", "Discordant differential expression",
                                  "No differential expression")))

all_log2fc_cor <- cor.test(c_to_d$B73_Log2FC, c_to_d$Mo17_Log2FC, method = "pearson")
significant_log2fc_cor <- cor.test(filter(c_to_d, B73_FDR < 0.05 & Mo17_FDR < 0.05)$B73_Log2FC, filter(c_to_d, B73_FDR < 0.05 & Mo17_FDR < 0.05)$Mo17_Log2FC, method = "pearson")

### Graph Stats ###
#################################################################################################################
scatter_plot_func(df = c_to_d, box_df = NULL, draw_boxes = FALSE, 
                  cor_value = round(all_log2fc_cor$estimate, digits = 2), 
                  p_value = "2.2e-16")

scatter_plot_func(df = filter(c_to_d, B73_FDR < 0.05 & Mo17_FDR < 0.05), draw_boxes = TRUE, 
                  box_df = boxes, cor_value = round(significant_log2fc_cor$estimate, digits = 2),
                  p_value = "2.2e-16")