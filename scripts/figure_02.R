### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(readxl)
library(DESeq2)
library(GO.db)
library(topGO)

### Load data ###
#################################################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
  mutate_at(vars(-Library_Name), factor)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name)

go_annotation_descriptions <- fread("./data/GO_Annotations.csv")

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% 
  distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% 
  mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

### Define functions ###
#################################################################################################################
# Differential expression analysis
DiffExp_func <- function(count_data, genotypes = c("B73", "Mo17"), sections = c("A", "B", "C", "D", "E"), silk_locations = c("Emerged", "Husk-Encased"), years = c("2014", "2015"), 
                         harvest_date = c("7/30/14", "7/30/15", "8/2/15", "7/31/15", "8/4/15"), design, contrast, read_depth_cutoff = 5, prop_libs_with_reads = 0.25, FDR = 0.05, FC = 2){
  
  sample_info <- filter(library_info, UQ(as.name(contrast[1])) %in% contrast[2:3], Genotype %in% genotypes, Section %in% sections, Year %in% years, Location_Of_Silk %in% silk_locations, 
                        Harvest_Date %in% Harvest_Date) %>%
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

# Gene enrichment of DE gene sets with TopGO
GO_enrichment_func <- function(myInterestingGenes, go_file_path, target_ontologies, p_value_cutoff){
  
  geneID2GO <- readMappings(file = go_file_path, sep = "\t", IDsep = ",")
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  ontologies <- target_ontologies
  for(j in seq_along(ontologies)){
    go_data <- new(Class = "topGOdata", ontology = ontologies[j], allGenes = geneList,
                   nodeSize = 10, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    Fisher_Test <- runTest(go_data, algorithm = "weight", statistic = "fisher")
    temp_go_results <- GenTable(go_data, P_Value = Fisher_Test, topNodes = length(Fisher_Test@score[Fisher_Test@score <= p_value_cutoff]))
    colnames(temp_go_results) <- c("GO", "Term", "Annotated", "Significant", "Expected", "P_Value")
    go_results <-  mutate(temp_go_results, P_Value = as.numeric(sub('<', '', P_Value)), 
                          Ontology_Type = as.character(ontologies[j]))
  }
  return(go_results)
}

### Compute Stats ###
#################################################################################################################
# Differential expression analysis across silk section transitions
sections <- c("A", "B", "C", "D", "E")
degs_across_sections <- ""
'%ni%' <- Negate('%in%')
for(i in 1:(length(sections) - 1)){
  b73_temp_df <- DiffExp_func(count_data = b73_counts, genotypes = "B73", design = ~ Section + Year, contrast = c("Section", sections[i], sections[i + 1])) %>% 
    mutate(Transition = paste(sections[i], "to", sections[i + 1], sep = " "))
  mo17_temp_df <- DiffExp_func(count_data = mo17_counts, genotypes = "Mo17", design = ~ Section + Year, contrast = c("Section", sections[i], sections[i + 1])) %>% 
    mutate(Transition = paste(sections[i], "to", sections[i + 1], sep = " "))
  
  b73_df <- left_join(dplyr::rename(b73_temp_df, "B73_V4_GeneID" = GeneID), syntenic_gene_pairs, by = "B73_V4_GeneID") %>% 
    mutate(Gene_Type = ifelse(Mo17_CAU_GeneID %in% mo17_temp_df$GeneID, "Commonly Differentially Expressed Syntelog", 
                              ifelse(Mo17_CAU_GeneID %ni% mo17_temp_df$GeneID & Mo17_CAU_GeneID %ni% syntenic_gene_pairs$Mo17_CAU_GeneID, "Non-syntenic", "Syntenic, genotype specific differential expression"))) %>%
    dplyr::select(Syntenic_Pair, B73_V4_GeneID, Mo17_CAU_GeneID, Transition, Expression_Change, Gene_Type)
  
  mo17_df <- left_join(dplyr::rename(mo17_temp_df, "Mo17_CAU_GeneID" = GeneID), syntenic_gene_pairs, by = "Mo17_CAU_GeneID") %>% 
    mutate(Gene_Type = ifelse(B73_V4_GeneID %in% b73_temp_df$GeneID, "Commonly Differentially Expressed Syntelog", 
                              ifelse(B73_V4_GeneID %ni% b73_temp_df$GeneID & B73_V4_GeneID %ni% syntenic_gene_pairs$B73_V4_GeneID, "Non-syntenic", "Syntenic, genotype specific differential expression"))) %>%
    dplyr::select(Syntenic_Pair, B73_V4_GeneID, Mo17_CAU_GeneID, Transition, Expression_Change, Gene_Type)
  
  syntenic_compare <- full_join(filter(b73_df, Gene_Type == "Commonly Differentially Expressed Syntelog") %>% dplyr::rename("B73_Expression_Change" = Expression_Change),
                                filter(mo17_df, Gene_Type == "Commonly Differentially Expressed Syntelog") %>% dplyr::rename("Mo17_Expression_Change" = Expression_Change), 
                                by = c("Syntenic_Pair", "B73_V4_GeneID", "Mo17_CAU_GeneID", "Gene_Type", "Transition")) %>%
    mutate(Gene_Type = ifelse(B73_Expression_Change == Mo17_Expression_Change, "Syntenic, concordant differential expression", "Syntenic, discordant differential expression"))
  
  b73_df <- bind_rows((dplyr::rename(b73_df, "GeneID" = B73_V4_GeneID) %>% mutate(Genotype = "B73") %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type) %>% filter(Gene_Type != "Commonly Differentially Expressed Syntelog")), 
                      (mutate(syntenic_compare, Genotype = "B73") %>% dplyr::rename("GeneID" = B73_V4_GeneID, "Expression_Change" = B73_Expression_Change) %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type)))
  
  mo17_df <- bind_rows((dplyr::rename(mo17_df, "GeneID" = Mo17_CAU_GeneID) %>% mutate(Genotype = "Mo17") %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type) %>% filter(Gene_Type != "Commonly Differentially Expressed Syntelog")), 
                       (mutate(syntenic_compare, Genotype = "Mo17") %>% dplyr::rename("GeneID" = Mo17_CAU_GeneID, "Expression_Change" = Mo17_Expression_Change) %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type)))
  
  degs_across_sections <- rbind(degs_across_sections, b73_df, mo17_df)
}
degs_across_sections_summary <- degs_across_sections[-1, ] %>% group_by(Genotype, Transition, Expression_Change, Gene_Type) %>% summarise(Gene_Count = n()) %>% 
  ungroup() %>% tidyr::complete(Genotype, Transition, Gene_Type, Expression_Change, fill = list(Gene_Count = 0)) %>%
  mutate(Gene_Count = ifelse(Expression_Change == "Expression Increases", Gene_Count, -1*Gene_Count))

# GO enrichment of syntenic C vs D DEGs with TopGO - biological processes; cellular component; molecular function
c_to_d_go_results <- tibble(GO = as.character(), Term = as.character(), 
                            Annotated = as.integer(), Significant = as.integer(), 
                            Expected = as.numeric(), P_Value = as.numeric(), 
                            Ontology_Type = as.character())
domains <- c("BP", "CC", "MF")
for(i in seq_along(domains)){
  c_to_d_go_results_temp <- GO_enrichment_func(myInterestingGenes = filter(degs_across_sections, Genotype == "B73", Transition == "C to D", Gene_Type == "Syntenic, concordant differential expression")$GeneID, 
                                               go_file_path = "./data/V4_B73_maize_GAMER_GOs.txt", 
                                               target_ontologies = domains[i], p_value_cutoff = 0.001)
  
  c_to_d_go_results_temp <- dplyr::select(c_to_d_go_results_temp, -Term) %>% left_join(., go_annotation_descriptions, by = "GO")
  
  c_to_d_go_results <- bind_rows(c_to_d_go_results, c_to_d_go_results_temp)
}

c_to_d_go_results <- c_to_d_go_results %>% group_by(Ontology_Type) %>% top_n(-20, P_Value) %>%
  mutate(Term = ifelse(Term == "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", 
                       "oxidoreductase activity, acting on paired donors", Term))
### Graph ###
#################################################################################################################
# Differential expression across section transitions
ggplot(filter(degs_across_sections_summary, Gene_Type != "Syntenic, discordant differential expression"), aes(x = Transition, y = Gene_Count, fill = reorder(Gene_Type, Gene_Count))) +
  guides(fill = FALSE) +
  geom_bar(stat = "identity", color = "black", size = 0.75, width = 0.75) +
  facet_grid(. ~ Genotype) +
  scale_fill_brewer(palette = "Greys") +
  geom_hline(yintercept = 0, size = 1.25, color = "black") +
  scale_y_continuous(limits = c(-500, 1750), breaks = c(-500, 0, 500, 1000, 1500), 
                     labels = c("500", "0", "500", "1000", "1500")) +
  theme_bw() +
  theme(plot.title = element_blank(), 
        axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), 
        panel.border = element_rect(size = 5),
        legend.text = element_text(size = 16), 
        legend.background = element_rect(color = "black", size = 1), 
        legend.position = "right", 
        strip.background = element_rect(fill = "white", size = 2, color = "black"), 
        strip.text = element_text(size = 64, face = "bold"), 
        panel.spacing = unit(0.05, "lines")) +
  labs(y = "# of Differentially\nExpressed Genes", fill = NULL)

# GO enrichment syntenic C vs D DEGs with TopGO
ggplot(c_to_d_go_results, aes(x = reorder(Term, Significant/Annotated), y = 100*Significant/Annotated, size = Annotated, color = -log10(P_Value))) +
  geom_point() + coord_flip() + facet_grid(Ontology_Type ~ ., scales = "free", space = "free_y") +
  scale_color_gradient(low = "blue", high = "red", breaks = c(4, 12, 20, 28), labels = c("4", "12", "20", "28")) +
  scale_size_continuous(breaks = c(100, 1000, 10000), range = c(2, 15)) +
  scale_y_continuous(limits = c(0, 31), breaks = c(0, 10, 20, 30), labels = c("0", "10", "20", "30")) +
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 24), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        legend.position = c(0.65, 0.7),
        legend.background = element_rect(size = 0.5, color = "black"),
        panel.border = element_rect(size = 3),
        strip.background = element_rect(fill = "white", size = 2, color = "black"), 
        strip.text = element_text(size = 40, face = "bold"),
        plot.margin = unit(c(1,1,1,1.25), "cm")) +
  labs(x = NULL,
       y = "% of Annotated Genes\nFound Significant", 
       size = "Annotated Genes\nper GO Term", 
       color = "-log10(p-value)")
