### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(DESeq2)

### Load data ###
#################################################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
                    mutate(Location_Of_Silk_Year = paste(Location_Of_Silk, "-", Year, sep = ""), 
                           group = paste(Section, Year, sep = ":")) %>%
                    mutate_at(vars(-Library_Name), factor)

b73_log2_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv")  %>%
                      dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name) %>%
                      mutate_at(vars(-GeneID), .funs = funs(log2((. + 1))))

mo17_log2_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
                      dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name) %>%
                      mutate_at(vars(-GeneID), .funs = funs(log2((. + 1))))

### Compute stats ###
#################################################################################################################
# Principal component analysis
pca_func <- function(log2_rpkm, genotype){
                log2_rpkm <- t(log2_rpkm %>% column_to_rownames(var = "GeneID"))
                sample_order <- rownames(log2_rpkm)
                
                syntenic_pca <- prcomp(log2_rpkm, center = TRUE)
                
                syntenic_pca_results <- tibble(PC1 = syntenic_pca$x[, 1], 
                                               PC2 = syntenic_pca$x[, 2], 
                                               PC3 = syntenic_pca$x[, 3], 
                                               Library_Name = sample_order) %>%
                                          left_join(., library_info, by = "Library_Name") %>%
                                          mutate(Genotype_Year = factor(paste(Genotype, ".", Year, sep = "")))
                
                percentVars <- syntenic_pca$sdev^2/sum(syntenic_pca$sdev^2)
                
                return(list("syntenic_pca_results" = syntenic_pca_results, "percentVars" = percentVars))
}

# PCA of B73 and Mo17 libraries separately
b73_pca_data <- pca_func(log2_rpkm = b73_log2_rpkm, genotype = "B73")
mo17_pca_data <- pca_func(log2_rpkm = mo17_log2_rpkm, genotype = "Mo17")

### Graph ###
#################################################################################################################
# PCA graphs
pca_plot_func <- function(pca_df, plot_title, annotation_label, percent_vars, x_pos, y_pos){
                    print(ggplot(pca_df, aes(PC1, PC2, color = Section, shape = Year)) +
                            geom_point(size = 8) +
                            scale_shape_manual(values = c(15, 17)) +
                            scale_x_continuous(limits = c(-134, 80), breaks = c(-125, -100, -75, -50, -25, 0, 25, 50, 75)) +
                            scale_y_continuous(limits = c(-80, 80), breaks = c(-75, -50, -25, 0, 25, 50, 75)) +
                            guides(shape = guide_legend(override.aes = list(size = 10)),
                                   color = guide_legend(override.aes = list(size = 10))) +
                            stat_ellipse(aes(PC1, PC2, group = Location_Of_Silk_Year), color = "black") +
                            xlab(paste0("PC1 (", round(100 * percent_vars[1], 0), "%)", sep = "")) +
                            ylab(paste0("PC2 (", round(100 * percent_vars[2], 0), "%)", sep = "")) +
                            coord_fixed() +
                            theme_bw() +
                            theme(plot.title = element_text(hjust = 0.5, size = 36), 
                                  axis.title = element_text(size = 32), 
                                  axis.text = element_text(size = 28), 
                                  legend.title = element_text(size = 32), 
                                  legend.text = element_text(size = 20), 
                                  panel.border = element_rect(size = 6)) +
                            ggtitle(plot_title) +
                            annotate(geom = "text", x = x_pos, y = y_pos, label = annotation_label, size = 12) +
                            labs(shape = "Year", color = "Section"))
}

pca_plot_func(pca_df = b73_pca_data$syntenic_pca_results, plot_title = "B73 Samples",
              annotation_label = paste("n = ", format(nrow(b73_log2_rpkm), big.mark = ",", scientific = FALSE), " genes"), 
              x_pos = -70, y_pos = 65, percent_vars = b73_pca_data$percentVars)

pca_plot_func(pca_df = mo17_pca_data$syntenic_pca_results, plot_title = "Mo17 Samples", 
              annotation_label = paste("n = ", format(nrow(mo17_log2_rpkm), big.mark = ",", scientific = FALSE), " genes"), 
              x_pos = -70, y_pos = 65, percent_vars = mo17_pca_data$percentVars)