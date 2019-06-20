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

tfs <- fread("./data/GRASSIUS_TFs.csv")

b73_expressed_genes <- fread("./data/RPKM_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, filter(library_info, Genotype == "B73")$Library_Name) %>%
  gather("Library_Name", "FPKM", 2:41) %>% 
  left_join(., dplyr::select(library_info, Library_Name, Section), by = "Library_Name") %>%
  group_by(GeneID, Section) %>% 
  summarise(FPKM = mean(FPKM, na.rm = TRUE)) %>% 
  filter(FPKM > 1)

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

tf_enrichment_func <- function(gene_list, p_value = 0.01){
  tf_deg_list <- list()
  tf_families <- unique(tfs$TF_Family)
  for(i in 1:length(tf_families)){
    temp_list <- list(filter(tfs, TF_Family == tf_families[i])$GeneID)
    tf_deg_list[tf_families[i]] <- temp_list
  }
  tf_deg_list["Differentially Expressed"] <- list(gene_list)
  tf.gom.obj <- newGOM(tf_deg_list, genome.size = 39320)
  tf_deg_overlap_pvals <- getMatrix(tf.gom.obj, "pval")
  tf_deg_overlap_pvals[lower.tri(tf_deg_overlap_pvals, diag = FALSE)] <- NA
  tf_deg_overlap_pvals <- melt(tf_deg_overlap_pvals, value.name = "p-value") %>% drop_na() %>%
    filter(!Var1 %in% c("Differentially Expressed"), 
           Var2 %in% c("Differentially Expressed"), 
           `p-value` < p_value)
  
  return(tf_deg_overlap_pvals)
}

### Compute Stats ###
#################################################################################################################
### Husk-encased vs emerged DEGs ###
encased_vs_emerged_degs <- DiffExp_func(b73_counts, "B73", design = ~ Location_Of_Silk + Year, contrast = c("Location_Of_Silk", "Husk-Encased", "Emerged")) %>%
  dplyr::select(GeneID, Expression_Change)

### TF enrichment tests ###
top_20_degs <- left_join(tfs, encased_vs_emerged_degs, by = "GeneID") %>%
  drop_na() %>%
  group_by(TF_Family) %>%
  summarise(Count = n()) %>%
  top_n(20, Count)

expressed_tfs <- filter(tfs, GeneID %in% unique(c(b73_expressed_genes$GeneID, encased_vs_emerged_degs$GeneID))) %>%
  group_by(TF_Family) %>%
  summarise(Expressed = n()) %>%
  left_join(., group_by(tfs, TF_Family) %>% summarise(Total = n()), by = "TF_Family") %>%
  mutate(`Not Expressed` = Total-Expressed) %>%
  dplyr::select(-Total) %>%
  gather("Gene_Type", "Count", 2:ncol(.)) %>%
  filter(TF_Family %in% top_20_degs$TF_Family) %>%
  mutate(TF_Family = reorder(TF_Family, Count), 
         Gene_Type = factor(Gene_Type, levels = c("Not Expressed", "Expressed")))

deg_tfs <- left_join(tfs, encased_vs_emerged_degs, by = "GeneID") %>%
  drop_na() %>%
  mutate(Gene_Type = ifelse(Expression_Change == "Expression Increases", "Up-regulated", "Down-regulated")) %>%
  group_by(TF_Family, Gene_Type) %>%
  summarise(Count = n()) %>%
  mutate(Gene_Type = factor(Gene_Type, levels = c("Up-regulated", "Down-regulated"))) %>%
  ungroup() %>%
  filter(TF_Family %in% top_20_degs$TF_Family)%>%
  mutate(TF_Family = factor(TF_Family, levels = levels(expressed_tfs$TF_Family)))

encased_vs_emerged_enriched_tfs <- tf_enrichment_func(gene_list = encased_vs_emerged_degs$GeneID, p_value = 0.01)

# Determine if enriched TF families in emerged silks are up- or down-regulated in seedling stress treatments
tf_stress_deg_direction <- left_join(stress_degs, tfs, by = "GeneID") %>%
  mutate(Expression_Change = ifelse(Expression_Change == "Decreases under treatment", "Down-regulated", "Up-regulated")) %>%
  group_by(Experiment, TF_Family, Expression_Change) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  filter(TF_Family %in% c("AP2-EREBP", "MYB", "NAC", "WRKY", "ZIM", "HSF")) %>%
  mutate(Expression_Change = factor(Expression_Change, levels = c("Up-regulated", "Down-regulated")), 
         TF_Family = factor(TF_Family, levels = c("HSF", "ZIM", "NAC", "WRKY", "MYB", "AP2-EREBP"))) %>%
  mutate(Experiment = ifelse(Experiment == "Cold stress", "Cold\nstress", 
                             ifelse(Experiment == "Heat stress", "Heat\nstress", 
                                    ifelse(Experiment == "Drought stress", "Drought\nstress", 
                                           ifelse(Experiment == "Salt stress", "Salt\nstress", 
                                                  ifelse(Experiment == "UV stress", "UV\nstress", 
                                                         ifelse(Experiment == "Jasmonic acid", "Jasmonic\nacid", 
                                                                ifelse(Experiment == "Salicylic acid", "Salicylic\nacid", "NA"))))))))

tf_stress_deg_direction[is.na(tf_stress_deg_direction)] <- 0

### Graph Stats ###
#################################################################################################################
ggplot(expressed_tfs, aes(x = TF_Family, y = -Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", color = "black", size = 1) + 
  coord_flip() +
  scale_y_continuous(breaks = c(0, -100, -200, -300), labels = c("0", "100", "200", "300"), expand = expand_scale(mult = c(0.05, 0))) +
  theme_bw() +
  scale_fill_manual(values = c("white", "grey15")) +
  theme(axis.text.x = element_text(size = 32),
        axis.title.x = element_text(size = 36),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 32), 
        legend.position = c(0.4, 0.25), 
        legend.background = element_rect(size = 1, color = "black"), 
        panel.border = element_rect(size = 3)) +
  labs(y = "Number of Genes", fill = NULL)

ggplot(deg_tfs, aes(x = TF_Family, y = Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", color = "black", size = 1) + 
  scale_y_continuous(limits = c(0, 60), expand = expand_scale(mult = c(0, .05))) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("firebrick3", "blue")) +
  theme(axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32, hjust = 0.5),
        axis.title.x = element_text(size = 36),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 32), 
        legend.position = c(0.6, 0.25), 
        legend.background = element_rect(size = 1, color = "black"),
        panel.border = element_rect(size = 4)) +
  labs(y = "Number of DEGs", fill = NULL)

ggplot(tf_stress_deg_direction, aes(x = TF_Family, y = Count, fill = Expression_Change)) +
  facet_grid(. ~ Experiment) +
  geom_bar(stat = "identity", color = "black", size = 0.75, width = 0.5) + 
  scale_y_continuous(limits = c(0, 40), expand = expand_scale(mult = c(0, .05))) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("firebrick3", "blue")) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 28),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 20), 
        legend.position = "bottom", 
        legend.background = element_rect(size = 1, color = "black"),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 28, 
                                  margin = margin(t = 0, r = 0, b = 10, l = 0)),
        panel.border = element_rect(size = 4), 
        panel.spacing = unit(0.75, "lines")) +
  labs(y = "Number of DEGs", fill = NULL)
