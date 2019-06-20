### Load libraries ###
#################################################################################
library(tidyverse)
library(readr)
library(data.table)
library(rlang)
library(WGCNA)
library(DESeq2)

### Load data and normalize ###
#################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
  mutate(Location_Of_Silk_Year = paste(Location_Of_Silk, "-", Year, sep = ""), 
         Genotype_Year = paste(Genotype, Year, sep = " "),
         group = paste(Section, Year, sep = ":")) %>% mutate_at(vars(-Library_Name), factor)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name)

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

combined_syntenic_module_assignments <- fread("./data/Combined_Syntenic_Network_Module_Assignments.csv") %>%
  left_join(., syntenic_gene_pairs, by = "Syntenic_Pair")

module_lookup_table <- tibble(Module = c('cyan', 'honeydew1', 'floralwhite', 'brown4', 'darkgrey', 'darkturquoise', 'darkolivegreen2', 
                                         'darkmagenta', 'darkred', 'darkseagreen3', 'midnightblue', 'darkviolet', 'firebrick4', 
                                         'lightsteelblue', 'yellowgreen', 'lightcyan1', 'white', 'mediumorchid', 'orangered4', 
                                         'greenyellow', 'darkslateblue', 'indianred4', 'firebrick3', 'honeydew', 'lavenderblush2', 'mediumpurple4', 'Unassigned'), 
                              New_Module = c(paste("Cluster", 1:26), "Unassigned"))

combined_syntenic_module_assignments <- left_join(combined_syntenic_module_assignments, module_lookup_table, by = "Module")

stress_degs <- fread("./data/Stress_DEGs.csv")

### Create functions ###
#################################################################################
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

### Compute stats ###
#################################################################################
### Differentially expressed gene analysis
encased_vs_emerged_degs <- bind_rows(DiffExp_func(b73_counts, "B73", design = ~ Location_Of_Silk + Year, contrast = c("Location_Of_Silk", "Husk-Encased", "Emerged")) %>%
                                      dplyr::select(GeneID, Expression_Change) %>% mutate(Genotype = "B73"),
                                    DiffExp_func(mo17_counts, "Mo17", design = ~ Location_Of_Silk + Year, contrast = c("Location_Of_Silk", "Husk-Encased", "Emerged")) %>%
                                      dplyr::select(GeneID, Expression_Change) %>% mutate(Genotype = "Mo17"))

### Total DEGs per cluster ###
degs <- bind_rows(left_join(filter(encased_vs_emerged_degs, Genotype == "B73"), syntenic_gene_pairs, by = c("GeneID" = "B73_V4_GeneID")) %>%
                    dplyr::select(Syntenic_Pair, Expression_Change, Genotype), 
                  left_join(filter(encased_vs_emerged_degs, Genotype == "Mo17"), syntenic_gene_pairs, by = c("GeneID" = "Mo17_CAU_GeneID")) %>%
                    dplyr::select(Syntenic_Pair, Expression_Change, Genotype)) %>%
  drop_na() %>% spread(Genotype, Expression_Change) %>%
  left_join(., combined_syntenic_module_assignments, by = "Syntenic_Pair") %>%
  filter(!is.na(New_Module)) %>%
  mutate(B73 = ifelse(is.na(B73), "NS", B73), 
         Mo17 = ifelse(is.na(Mo17), "NS", Mo17)) %>%
  mutate(Gene_Type_1 = ifelse(!B73 %in% c("Expression Increases", "Expression Decreases"), "Mo17 Only",
                              ifelse(!Mo17 %in% c("Expression Increases", "Expression Decreases"), "B73 Only", "B73 & Mo17")),
         Gene_Type_2 = ifelse(B73 == "Expression Increases" & Mo17 == "Expression Decreases", "Opposite", 
                              ifelse(B73 == "Expression Decreases" & Mo17 == "Expression Increases", "Opposite",
                                     ifelse(B73 == "Expression Increases" & Mo17 == "Expression Increases", "Up-regulated",
                                            ifelse(B73 == "Expression Decreases" & Mo17 == "Expression Decreases", "Down-regulated",
                                                   ifelse(B73 == "Expression Increases" & Mo17 == "NS", "Up-regulated",
                                                          ifelse(Mo17 == "Expression Increases" & B73 == "NS", "Up-regulated",
                                                                 ifelse(B73 == "Expression Decreases" & Mo17 == "NS", "Down-regulated",
                                                                        ifelse(Mo17 == "Expression Decreases" & B73 == "NS", "Down-regulated", "NA"))))))))) %>%
  group_by(New_Module, Gene_Type_1, Gene_Type_2) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(New_Module = factor(New_Module, levels = c("Cluster 10",
                                                    "Cluster 14",
                                                    "Cluster 16",
                                                    "Cluster 13",
                                                    "Cluster 15",
                                                    "Cluster 11",
                                                    "Cluster 17",
                                                    "Cluster 22",
                                                    "Cluster 18",
                                                    "Cluster 12",
                                                    "Cluster 24",
                                                    "Cluster 23",
                                                    "Cluster 26",
                                                    "Cluster 8",
                                                    "Cluster 2",
                                                    "Cluster 25",
                                                    "Cluster 21",
                                                    "Cluster 5",
                                                    "Cluster 9",
                                                    "Cluster 3",
                                                    "Cluster 1",
                                                    "Cluster 7",
                                                    "Cluster 20",
                                                    "Cluster 19",
                                                    "Cluster 6",
                                                    "Cluster 4")))

stress_deg_list <- list("Cold Up" = filter(stress_degs, Experiment == "Cold stress", Expression_Change == "Increases under treatment")$GeneID,
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

# unique syntenic genes up-regulated under any one of the stresses and also up-regulated in emerged silks
syn_upregulated_genes <- intersect(syntenic_gene_pairs$B73_V4_GeneID, intersect(unique(c(stress_deg_list$`Cold Up`, 
                                                                                         stress_deg_list$`Heat Up`, 
                                                                                         stress_deg_list$`Drought Up`, 
                                                                                         stress_deg_list$`Salt Up`, 
                                                                                         stress_deg_list$`UV Up`,
                                                                                         stress_deg_list$`JA Up`,
                                                                                         stress_deg_list$`SA Up`)), stress_deg_list$`Encased to Emerged Up`))

# overlap of up-regulated syntenic genes and genes in clusters 10, 14, and 16
length(intersect(syn_upregulated_genes, 
                 filter(combined_syntenic_module_assignments, 
                        New_Module %in% c("Cluster 10", "Cluster 14", "Cluster 16"))$B73_V4_GeneID))

# unique syntenic genes down-regulated in drought stress and also down-regulated in emerged silks
syn_downregulated_genes <- intersect(syntenic_gene_pairs$B73_V4_GeneID, intersect(stress_deg_list$`Drought Down`, 
                                                                                  stress_deg_list$`Encased to Emerged Down`))

# overlap of down-regulated syntenic genes and genes in clusters 4 and 6
length(intersect(syn_downregulated_genes, 
                 filter(combined_syntenic_module_assignments, 
                        New_Module %in% c("Cluster 4", "Cluster 6"))$B73_V4_GeneID))

### Graph summary stats of combined syntenic network ###
#################################################################################
### Differentially expressed gene counts ###
ggplot(filter(degs, Gene_Type_2 == "Down-regulated", New_Module != "Unassigned"), aes(x = New_Module, y = -Count, fill = Gene_Type_1)) + 
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(. ~ Gene_Type_2) +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-165, 0), breaks = c(-100, 0), labels = c(100, 0), expand = expand_scale(mult = c(0, 0))) +
  coord_flip() +
  theme(axis.title.x = element_text(size = 32), 
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 36, 
                                  margin = margin(t = 0, r = 0, b = 20, l = 0)), 
        panel.border = element_rect(size = 3),
        legend.position = "none") +
  labs(x = NULL, fill = NULL, y = "# of Genes")

ggplot(filter(degs, Gene_Type_2 == "Up-regulated", New_Module != "Unassigned"), aes(x = New_Module, y = Count, fill = Gene_Type_1)) + 
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(. ~ Gene_Type_2) +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 575), expand = expand_scale(mult = c(0, .05))) +
  coord_flip() +
  theme(axis.title = element_text(size = 32), 
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 36, 
                                  margin = margin(t = 0, r = 0, b = 20, l = 0)),
        plot.margin = unit(c(1, 1.25, 1, 1),"cm"),
        legend.title = element_text(size = 32), 
        legend.text = element_text(size = 28), 
        panel.border = element_rect(size = 3), 
        legend.position = c(0.6, 0.5), 
        legend.background = element_rect(size = 1.5, color = "black")) +
  labs(x = NULL, fill = "Differential expression\nbetween emerged\nand encased in:", y = "# of Genes")