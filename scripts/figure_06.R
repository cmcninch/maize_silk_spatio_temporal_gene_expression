### Load libraries ###
#################################################################################
library(tidyverse)
library(readr)
library(data.table)
library(rlang)
library(WGCNA)
library(DESeq2)
library(GO.db)
library(topGO)

### Load data and normalize ###
#################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
  mutate(Location_Of_Silk_Year = paste(Location_Of_Silk, "-", Year, sep = ""), 
         Genotype_Year = paste(Genotype, Year, sep = " "),
         group = paste(Section, Year, sep = ":")) %>% mutate_at(vars(-Library_Name), factor)

b73_expressed_genes <- fread("./data/RPKM_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, filter(library_info, Genotype == "B73")$Library_Name) %>%
  gather("Library_Name", "RPKM", 2:41) %>% left_join(., dplyr::select(library_info, Library_Name, Section), by = "Library_Name") %>%
  group_by(GeneID, Section) %>% summarise(RPKM = mean(RPKM, na.rm = TRUE)) %>% filter(RPKM > 1) %>% distinct(GeneID)

mo17_expressed_genes <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, filter(library_info, Genotype == "Mo17")$Library_Name) %>%
  gather("Library_Name", "RPKM", 2:41) %>% left_join(., dplyr::select(library_info, Library_Name, Section), by = "Library_Name") %>%
  group_by(GeneID, Section) %>% summarise(RPKM = mean(RPKM, na.rm = TRUE)) %>% filter(RPKM > 1) %>% distinct(GeneID)

b73_log2_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv")  %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name) %>%
  mutate_at(vars(-GeneID), .funs = funs(log2((. + 1)))) %>% filter(GeneID %in% b73_expressed_genes$GeneID)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name)

mo17_log2_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name) %>%
  mutate_at(vars(-GeneID), .funs = funs(log2((. + 1)))) %>% filter(GeneID %in% mo17_expressed_genes$GeneID)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name)

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

log2_rpkm <- left_join(filter(b73_log2_rpkm, GeneID %in% syntenic_gene_pairs$B73_V4_GeneID), syntenic_gene_pairs, by = c("GeneID" = "B73_V4_GeneID")) %>%
  dplyr::select(-GeneID) %>%
  left_join(., (left_join(filter(mo17_log2_rpkm, GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID), syntenic_gene_pairs, by = c("GeneID" = "Mo17_CAU_GeneID")) %>% dplyr::select(-GeneID)), 
            by = "Syntenic_Pair") %>%
  dplyr::select(-B73_V4_GeneID, -Mo17_CAU_GeneID) %>% dplyr::select(Syntenic_Pair, everything()) %>% 
  drop_na()

rownames(log2_rpkm) <- NULL

go_annotation_descriptions <- fread("./data/GO_Annotations.csv")

### Create functions ###
#################################################################################
sft_threshold_optimize_func <- function(rpkm_data){
  rpkm_data <- t(rpkm_data %>% column_to_rownames(var = "Syntenic_Pair"))
  sft = pickSoftThreshold(rpkm_data, powerVector = seq(2,30, by = 2), networkType = "signed", verbose = 0, corFnc = "bicor")
  # Graph the scale-free topology fit index and mean connectivity as a function of the soft-thresholding power
  sizeGrWindow(9, 5); par(mfrow = c(1,2)); cex1 = 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", 
       type = "n", main = paste("Scale independence"));
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = seq(2, 30, by = 2), cex = cex1, col = "red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", 
       ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = seq(2, 30, by = 2), 
       cex = cex1, col = "red")
  
  sft_threshold <- readline("What Soft Threshold Power is best? ")
  return(as.numeric(sft_threshold))
}

# Gene enrichment of cluster gene sets with TopGO
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

### Determine best soft threshold power ###
#################################################################################
sft_threshold <- sft_threshold_optimize_func(rpkm_data = log2_rpkm)

### Build network ###
#################################################################################
adjacency_matrix <- adjacency(datExpr = t(log2_rpkm %>% column_to_rownames(var = "Syntenic_Pair")),
                              type = "signed", power = 16, corFnc = "bicor",
                              corOptions = list(maxPOutliers = 0.1))

TOM_similarity_matrix <- TOMsimilarity(adjMat = adjacency_matrix)
TOM_distance_matrix = 1 - TOM_similarity_matrix

gene_tree = hclust(as.dist(TOM_distance_matrix), method = "average")

modules = cutreeDynamic(dendro = gene_tree, method = "tree", deepSplit = 2, minClusterSize = 30)

dynamic_colors = labels2colors(modules)

module_eigengene_list = moduleEigengenes(expr = t(log2_rpkm %>% column_to_rownames(var = "Syntenic_Pair")),
                                         colors = dynamic_colors)
module_eigengenes = module_eigengene_list$eigengenes

module_eigengenes_distance_matrix <-  1 - cor(module_eigengenes)

module_eigengenes_tree = hclust(as.dist(module_eigengenes_distance_matrix), method = "average")

merged_modules = mergeCloseModules(exprData = t(log2_rpkm %>% column_to_rownames(var = "Syntenic_Pair")),
                                   dynamic_colors, cutHeight = 0.15, verbose = 3)

merged_module_assignments <- tibble(Syntenic_Pair = log2_rpkm$Syntenic_Pair,
                                    Module = merged_modules$colors) %>% mutate(Module = ifelse(Module == "grey", "Unassigned", Module))
merged_module_eigengenes = merged_modules$newMEs %>% mutate("Library_Name" = c(filter(library_info, Genotype == "B73")$Library_Name,
                                                                               filter(library_info, Genotype == "Mo17")$Library_Name)) %>%
  gather("Module", "Eigengene_Expression_Level:", 1:(ncol(.) - 1)) %>%
  mutate(Module = gsub("ME", "", Module)) %>%
  left_join(., dplyr::select(library_info, Library_Name, Section), by = "Library_Name") %>%
  mutate(Module = ifelse(Module == "grey", "Unassigned", Module))

merged_module_connectivities <- intramodularConnectivity(adjacency_matrix, merged_modules$colors) %>%
  rownames_to_column(var = "Syntenic_Pair")

merged_module_assignments <- left_join(merged_module_assignments, merged_module_connectivities, by = "Syntenic_Pair")

### Obtain adjacency scores within each module ###
adjacency_scores <- reshape2::melt(replace(adjacency_matrix, lower.tri(adjacency_matrix, diag = TRUE), NA),
                                   na.rm = TRUE, varnames = c("Node1", "Node2"), value.name = "Adjacency_Score") %>%
  left_join(., merged_module_assignments, by = c("Node1" = "Syntenic_Pair")) %>%
  left_join(., merged_module_assignments, by = c("Node2" = "Syntenic_Pair")) %>%
  dplyr::rename("Node1_Module" = "Module.x", "Node2_Module" = "Module.y") %>%
  filter(Node1_Module == Node2_Module) %>%
  dplyr::select(Node1, Node2, Adjacency_Score, Node1_Module, Node2_Module)

### Save results ###
#################################################################################
write.csv(merged_module_assignments, "./data/Combined_Syntenic_Network_Module_Assignments.csv", row.names = FALSE, quote = FALSE)
write.csv(merged_module_connectivities, "./data/Combined_Syntenic_Network_Intramodular_Connectivity_Scores.csv", row.names = FALSE, quote = FALSE)
write.csv(adjacency_scores, "./data/Combined_Syntenic_Network_Adjacency_Scores.csv", row.names = FALSE, quote = FALSE)
write.csv(merged_module_eigengenes, "./data/Combined_Syntenic_Network_Module_Eigengenes.csv", row.names = FALSE, quote = FALSE)

### Load results ###
#################################################################################
combined_syntenic_module_assignments <- fread("./data/Combined_Syntenic_Network_Module_Assignments.csv") %>%
  left_join(., syntenic_gene_pairs, by = "Syntenic_Pair")
combined_syntenic_module_eigengenes <- fread("./data/Combined_Syntenic_Network_Module_Eigengenes.csv")

### Compute stats ###
#################################################################################
### Relationship amongst modules ###
merged_module_dist <- 1 - cor(dplyr::select(combined_syntenic_module_eigengenes, -Section) %>%
                                spread(Module, `Eigengene_Expression_Level:`) %>%
                                column_to_rownames(var = "Library_Name"))

merged_module_tree <- hclust(as.dist(merged_module_dist), method = "average")

### Eigengene averages ###
module_eigengene_averages <- left_join(combined_syntenic_module_eigengenes, dplyr::select(library_info, - Section), by = "Library_Name") %>%
  group_by(Module, Genotype, Section) %>% summarise(Expression_Mean = mean(`Eigengene_Expression_Level:`, na.rm = TRUE)) %>%
  ungroup()

module_lookup_table <- tibble(Module = c('cyan', 'honeydew1', 'floralwhite', 'brown4', 'darkgrey', 'darkturquoise', 'darkolivegreen2', 
                                         'darkmagenta', 'darkred', 'darkseagreen3', 'midnightblue', 'darkviolet', 'firebrick4', 
                                         'lightsteelblue', 'yellowgreen', 'lightcyan1', 'white', 'mediumorchid', 'orangered4', 
                                         'greenyellow', 'darkslateblue', 'indianred4', 'firebrick3', 'honeydew', 'lavenderblush2', 'mediumpurple4', 'Unassigned'), 
                              New_Module = c(paste("Cluster", 1:26), "Unassigned"))

module_eigengene_averages <- left_join(module_eigengene_averages, module_lookup_table, by = "Module") %>%
  mutate(New_Module = factor(New_Module, levels = c(paste("Cluster", 1:26), "Unassigned")))

### Gene counts per module ###
module_counts <- group_by(combined_syntenic_module_assignments, Module) %>%
  summarise(Total = n()) %>% 
  left_join(., module_lookup_table, by = "Module") %>%
  mutate(New_Module = factor(New_Module, levels = c(paste("Cluster", 1:26), "Unassigned")))

### GO enrichment of syntenic C vs D DEGs with TopGO - biological processes; cellular component; molecular function ###
combined_syntenic_module_assignments <- left_join(combined_syntenic_module_assignments, module_lookup_table, by = "Module")


cluster_go_results <- tibble(GO = as.character(), Term = as.character(), 
                             Annotated = as.integer(), Significant = as.integer(), 
                             Expected = as.numeric(), P_Value = as.numeric(), 
                             Ontology_Type = as.character(), Cluster = as.character())
clusters <- unique(combined_syntenic_module_assignments$New_Module)

for(i in seq_along(clusters)){
  domains <- c("BP", "CC", "MF")
  for(j in seq_along(domains)){
    go_results_temp <- GO_enrichment_func(myInterestingGenes = filter(combined_syntenic_module_assignments, New_Module == clusters[i])$B73_V4_GeneID, 
                                          go_file_path = "./data/V4_B73_maize_GAMER_GOs_Syntenic_Only.txt", 
                                          target_ontologies = domains[j], p_value_cutoff = 0.001) %>%
      mutate(Cluster = clusters[i])
    
    go_results_temp <- dplyr::select(go_results_temp, -Term) %>% left_join(., go_annotation_descriptions, by = "GO")
    
    cluster_go_results <- bind_rows(cluster_go_results, go_results_temp)
  }
}

### Graph summary stats of combined syntenic network ###
#################################################################################
### Module eigengene averages ###
ggplot(module_eigengene_averages, aes(x = Section, y = Expression_Mean, color = Genotype)) + 
  facet_wrap(vars(New_Module), nrow = 3) + geom_point(size = 3) + geom_line(aes(group = Genotype), size = 1.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 32), 
        axis.title.y = element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = 24),
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 20), 
        panel.border = element_rect(size = 3), 
        legend.position = "bottom", 
        legend.background = element_rect(size = 1.5, color = "black"), 
        strip.background = element_rect(fill = "orange", size = 2), 
        strip.text = element_text(size = 16)) +
  labs(x = "Section", fill = "Genotype", y = "Eigengene Mean Expression") + 
  geom_text(data = dplyr::select(module_counts, New_Module, Total), aes(x = 2.5, y = 0.2, label = paste("n = ", Total, sep = "")), 
            colour = "black", size = 6, inherit.aes = FALSE, parse = FALSE)