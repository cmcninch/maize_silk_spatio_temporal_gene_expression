### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(DESeq2)
library(GeneOverlap)
library(FactoMineR)
library(ggrepel)

### Load data ###
#################################################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
  mutate_at(vars(-Library_Name), factor)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name)

stress_degs <- fread("./data/Stress_DEGs.csv")

### Define functions ###
#################################################################################################################
# Differential expression analysis
DiffExp_func <- function(count_data, genotypes = c("B73", "Mo17"), sections = c("A", "B", "C", "D", "E"), silk_locations = c("Emerged", "Husk-Encased"), years = c("2014", "2015"), 
                         design, contrast, read_depth_cutoff = 5, prop_libs_with_reads = 0.25, FDR = 0.05, FC = 2){
  
  sample_info <- filter(library_info, UQ(as.name(contrast[1])) %in% contrast[2:3], Genotype %in% genotypes, Section %in% sections, Year %in% years, Location_Of_Silk %in% silk_locations) %>%
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

### Compute Stats ###
#################################################################################################################
### Stress DEG comparisons ###
encased_vs_emerged_degs <- DiffExp_func(b73_counts, "B73", design = ~ Location_Of_Silk + Year, contrast = c("Location_Of_Silk", "Husk-Encased", "Emerged")) %>%
  dplyr::select(GeneID, Expression_Change)

deg_list <- list("Cold Up" = filter(stress_degs, Experiment == "Cold stress", Expression_Change == "Increases under treatment")$GeneID,
                 "Cold Down" = filter(stress_degs, Experiment == "Cold stress", Expression_Change == "Decreases under treatment")$GeneID,
                 "Heat Up" = filter(stress_degs, Experiment == "Heat stress", Expression_Change == "Increases under treatment")$GeneID,
                 "Heat Down" = filter(stress_degs, Experiment == "Heat stress", Expression_Change == "Decreases under treatment")$GeneID,
                 "Drought Up" = filter(stress_degs, Experiment == "Drought stress", Expression_Change == "Increases under treatment")$GeneID,
                 "Drought Down" = filter(stress_degs, Experiment == "Drought stress", Expression_Change == "Decreases under treatment")$GeneID,
                 "Salt Up" = filter(stress_degs, Experiment == "Salt stress", Expression_Change == "Increases under treatment")$GeneID,
                 "Salt Down" = filter(stress_degs, Experiment == "Salt stress", Expression_Change == "Decreases under treatment")$GeneID,
                 "UV Up" = filter(stress_degs, Experiment == "UV stress", Expression_Change == "Increases under treatment")$GeneID,
                 "UV Down" = filter(stress_degs, Experiment == "UV stress", Expression_Change == "Decreases under treatment")$GeneID,
                 "JA Up" = filter(stress_degs, Experiment == "Jasmonic acid", Expression_Change == "Increases under treatment")$GeneID, 
                 "JA Down" = filter(stress_degs, Experiment == "Jasmonic acid", Expression_Change == "Decreases under treatment")$GeneID,
                 "SA Up" = filter(stress_degs, Experiment == "Salicylic acid", Expression_Change == "Increases under treatment")$GeneID,
                 "SA Down" = filter(stress_degs, Experiment == "Salicylic acid", Expression_Change == "Decreases under treatment")$GeneID,
                 "Encased to Emerged Up" = filter(encased_vs_emerged_degs, Expression_Change == "Expression Increases")$GeneID,
                 "Encased to Emerged Down" = filter(encased_vs_emerged_degs, Expression_Change == "Expression Decreases")$GeneID)

gom.obj <- newGOM(deg_list, genome.size = 39320)

deg_overlap_counts <- getMatrix(gom.obj, "intersection")
deg_overlap_counts[lower.tri(deg_overlap_counts, diag = FALSE)] <- NA
deg_overlap_counts <- melt(deg_overlap_counts, value.name = "Overlap_Size") %>% drop_na() %>%
  filter(!Var1 %in% c("Encased to Emerged Up", "Encased to Emerged Down"), 
         Var2 %in% c("Encased to Emerged Up", "Encased to Emerged Down"))

deg_overlap_p_values <- getMatrix(gom.obj, "pval")
deg_overlap_p_values[lower.tri(deg_overlap_p_values, diag = FALSE)] <- NA
deg_overlap_p_values <- melt(deg_overlap_p_values, value.name = "-log10(p-value)") %>% drop_na() %>%
  filter(!Var1 %in% c("Encased to Emerged Up", "Encased to Emerged Down"), 
         Var2 %in% c("Encased to Emerged Up", "Encased to Emerged Down")) %>%
  mutate(`-log10(p-value)` = round(-log10(`-log10(p-value)`), digits = 2))

deg_overlap_summary <- left_join(deg_overlap_counts, deg_overlap_p_values, by = c("Var1", "Var2")) %>%
  separate(Var1, c("Treatment", "Expression_Change")) %>%
  mutate(Expression_Change = ifelse(Expression_Change == "Up", "Up-regulated", "Down-regulated")) %>%
  mutate(Var2 = ifelse(Var2 == "Encased to Emerged Up", "Up", "Down")) %>%
  mutate(Treatment = ifelse(Treatment == "Cold", "Cold\nstress", 
                            ifelse(Treatment == "Heat", "Heat\nstress", 
                                   ifelse(Treatment == "Drought", "Drought\nstress", 
                                          ifelse(Treatment == "Salt", "Salt\nstress", 
                                                 ifelse(Treatment == "UV", "UV\nstress", 
                                                        ifelse(Treatment == "JA", "Jasmonic\nacid", 
                                                               ifelse(Treatment == "SA", "Salicylic\nacid", "NA")))))))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Drought\nstress", "Jasmonic\nacid", "Heat\nstress","Salicylic\nacid",
                                                  "UV\nstress", "Salt\nstress", "Cold\nstress")))

# number of unique genes up-regulated under any one of the stresses and also up-regulated in emerged silks
length(intersect(unique(c(deg_list$`Cold Up`, 
                          deg_list$`Heat Up`, 
                          deg_list$`Drought Up`, 
                          deg_list$`Salt Up`, 
                          deg_list$`UV Up`,
                          deg_list$`JA Up`,
                          deg_list$`SA Up`)), deg_list$`Encased to Emerged Up`))

# MCA of stresses and silks
stress_degs_spread <- mutate(stress_degs, Expression_Change = gsub(" under treatment", "", Expression_Change)) %>%
  distinct() %>% spread(Experiment, Expression_Change) %>%
  full_join(., encased_vs_emerged_degs, by = "GeneID") %>%
  mutate(`Silk Emergence` = gsub("Expression ", "", Expression_Change)) %>%
  dplyr::select(-Expression_Change) %>%
  distinct()

mca_results <- MCA(stress_degs_spread %>% column_to_rownames(var = "GeneID"), graph = FALSE)
percent_vars <- as.data.frame(mca_results$eig)$`percentage of variance`

res.hcpc <- HCPC(as.data.frame(mca_results$var$eta2), nb.clust = -1, graph = FALSE)$data.clust %>%
  rownames_to_column(var = "Treatment") %>%
  dplyr::rename(Cluster = clust) %>%
  dplyr::select(Treatment, Cluster)

mca_results_final <- as.data.frame(mca_results$var$eta2) %>% 
  rownames_to_column(var = "Treatment") %>%
  left_join(., res.hcpc, by = "Treatment") %>%
  mutate(Cluster = ifelse(Cluster == "1", "One", 
                          ifelse(Cluster == "2", "Two", "Three"))) %>%
  mutate(Cluster = factor(Cluster, levels = c("One", "Two", "Three")))

### Create figures ###
#################################################################################################################
ggplot(deg_overlap_summary, aes(x = Treatment, y = Overlap_Size, fill = `-log10(p-value)`, shape = Expression_Change)) +
  geom_point(size = 10) + 
  facet_grid(Var2 ~ .) + 
  scale_fill_continuous(low = "blue", high = "red",
                        labels = c("25", "75", "125", "175", "225"), 
                        breaks = c(25, 75, 125, 175, 225)) +
  guides(shape = guide_legend(override.aes = list(size = 6.5))) +
  scale_shape_manual(values = c(25, 24)) +
  theme_bw() +
  theme(axis.text = element_text(size = 24),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 30, 
                                  margin = margin(t = 0, r = 0, b = 20, l = 0)), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 16), 
        legend.key = element_rect(size = 16),
        legend.key.size = unit(1.75, 'lines'), 
        legend.position = "right", 
        legend.background = element_blank(), 
        panel.border = element_rect(size = 4), 
        legend.box = "vertical") +
  labs(shape = "Expression Change\nUnder Treatment", 
       fill = "-log10(p-value)", 
       x = NULL, 
       y = "# of Shared DEGs in Meta-analysis")

ggplot(mca_results_final, aes(x = `Dim 1`, y = `Dim 2`, color = Cluster)) +
  geom_point(size = 8) +
  scale_y_continuous(limits = c(-0.1, 0.65)) +
  geom_text_repel(aes(label = Treatment, color = Cluster), size = 9, 
                  point.padding = 0.5, segment.color = "transparent", show.legend = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 40, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 40, margin = margin (0, 10, 0, 0)),
        axis.text = element_text(size = 32), 
        panel.border = element_rect(size = 4), 
        legend.position = "none") +
  labs(x = paste0("Dim 1 (", round(percent_vars[1], 1),"% variance)"),
       y = paste0("Dim 2 (", round(percent_vars[2], 1),"% variance)"))