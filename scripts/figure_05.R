### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(readxl)
library(DESeq2)
library(GO.db)
library(topGO)
library(UpSetR)

### Load data ###
#################################################################################################################
library_codes <- fread("./data/Library_Codes.csv") %>%
  mutate(Location_Of_Silk_Year = paste(Location_Of_Silk, "-", Year, sep = ""), 
         group = paste(Section, Year, sep = ":")) %>%
  mutate_at(vars(-Library_Name), factor)

counts_with_b73_genome <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv")

counts_with_mo17_genome <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv")

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

go_annotation_descriptions <- fread("./data/GO_Annotations.csv")

### Define functions ###
#################################################################################################################
# Differential expression analysis
DiffExp_func <- function(count_data, genotypes = c("B73", "Mo17"), sections = c("A", "B", "C", "D", "E"), silk_locations = c("Emerged", "Husk-Encased"), years = c("2014", "2015"), 
                         design, contrast, read_depth_cutoff = 5, prop_libs_with_reads = 0.25, FDR = 0.05, FC = 2){
  
  sample_info <- filter(library_codes, UQ(as.name(contrast[1])) %in% contrast[2:3], Genotype %in% genotypes, Section %in% sections, Year %in% years, Location_Of_Silk %in% silk_locations) %>%
    arrange(UQ(as.name(contrast[1]))) %>% column_to_rownames(var = "Library_Name")
  count_data <- dplyr::select(count_data, GeneID, rownames(sample_info)) %>% arrange(GeneID)
  max_low_expressed_libraries <- (1 - prop_libs_with_reads) * (ncol(count_data) - 1)
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

### Differential expression analysis of each individual silk section between genotypes when aligned to and quantified with the B73 genome ###
#################################################################################################################
sections <- c("A", "B", "C", "D", "E")
b73_degs_between_genotypes <- tibble()

for(i in 1:(length(sections))){
  temp_results <- DiffExp_func(count_data = counts_with_b73_genome, sections = sections[i], design = ~ Genotype + Year, contrast = c("Genotype", "B73", "Mo17")) %>% 
    mutate(Section = sections[i])
  
  b73_degs_between_genotypes <- bind_rows(b73_degs_between_genotypes, temp_results)
}

### Differential expression analysis each individual silk section between genotypes when aligned to and quantified with the Mo17 genome ###
#################################################################################################################
sections <- c("A", "B", "C", "D", "E")
mo17_degs_between_genotypes <- tibble()

for(i in 1:(length(sections))){
  temp_results <- DiffExp_func(count_data = counts_with_mo17_genome, sections = sections[i], design = ~ Genotype + Year, contrast = c("Genotype", "B73", "Mo17")) %>% 
    mutate(Section = sections[i])
  
  mo17_degs_between_genotypes <- bind_rows(mo17_degs_between_genotypes, temp_results)
}

### Summarise how syntenic genes are differentially expressed across the two used genomes in each individual section###
#################################################################################################################
b73_deg_summary <- dplyr::select(b73_degs_between_genotypes, GeneID, Section, Expression_Change) %>%
  spread(Section, Expression_Change) %>% 
  left_join(filter(syntenic_gene_pairs, B73_V4_GeneID %in% b73_degs_between_genotypes$GeneID), ., by = c("B73_V4_GeneID" = "GeneID")) %>%
  dplyr::select(-Mo17_CAU_GeneID)
colnames(b73_deg_summary) <- c("B73_V4_GeneID", "Syntenic_Pair", "B73_A", "B73_B", "B73_C", "B73_D", "B73_E")

mo17_deg_summary <- dplyr::select(mo17_degs_between_genotypes, GeneID, Section, Expression_Change) %>%
  spread(Section, Expression_Change) %>% 
  left_join(filter(syntenic_gene_pairs, Mo17_CAU_GeneID %in% mo17_degs_between_genotypes$GeneID), ., by = c("Mo17_CAU_GeneID" = "GeneID")) %>%
  dplyr::select(-B73_V4_GeneID)
colnames(mo17_deg_summary) <- c("Mo17_CAU_GeneID", "Syntenic_Pair", "Mo17_A", "Mo17_B", "Mo17_C", "Mo17_D", "Mo17_E")

deg_summary <- full_join(b73_deg_summary, mo17_deg_summary, by = "Syntenic_Pair") %>%
  mutate(A = ifelse(B73_A == Mo17_A & !is.na(B73_A) & !is.na(Mo17_A), "Concordant", 
                    ifelse(is.na(B73_A) & is.na(Mo17_A), "No differential expression when using either genome", 
                           ifelse(is.na(B73_A) & !is.na(Mo17_A), "Mo17\ngenome\nonly", 
                                  ifelse(!is.na(B73_A) & is.na(Mo17_A), "B73\ngenome\nonly", 
                                         ifelse(B73_A != Mo17_A, "Discordant", "NA")))))) %>%
  mutate(B = ifelse(B73_B == Mo17_B & !is.na(B73_B) & !is.na(Mo17_B), "Concordant", 
                    ifelse(is.na(B73_B) & is.na(Mo17_B), "No differential expression when using either genome", 
                           ifelse(is.na(B73_B) & !is.na(Mo17_B), "Mo17\ngenome\nonly", 
                                  ifelse(!is.na(B73_B) & is.na(Mo17_B), "B73\ngenome\nonly", 
                                         ifelse(B73_B != Mo17_B, "Discordant", "NA")))))) %>%
  mutate(C = ifelse(B73_C == Mo17_C & !is.na(B73_C) & !is.na(Mo17_C), "Concordant", 
                    ifelse(is.na(B73_C) & is.na(Mo17_C), "No differential expression when using either genome", 
                           ifelse(is.na(B73_C) & !is.na(Mo17_C), "Mo17\ngenome\nonly", 
                                  ifelse(!is.na(B73_C) & is.na(Mo17_C), "B73\ngenome\nonly", 
                                         ifelse(B73_C != Mo17_C, "Discordant", "NA")))))) %>%
  mutate(D = ifelse(B73_D == Mo17_D & !is.na(B73_D) & !is.na(Mo17_D), "Concordant", 
                    ifelse(is.na(B73_D) & is.na(Mo17_D), "No differential expression when using either genome", 
                           ifelse(is.na(B73_D) & !is.na(Mo17_D), "Mo17\ngenome\nonly", 
                                  ifelse(!is.na(B73_D) & is.na(Mo17_D), "B73\ngenome\nonly", 
                                         ifelse(B73_D != Mo17_D, "Discordant", "NA")))))) %>%
  mutate(E = ifelse(B73_E == Mo17_E & !is.na(B73_E) & !is.na(Mo17_E), "Concordant", 
                    ifelse(is.na(B73_E) & is.na(Mo17_E), "No differential expression when using either genome", 
                           ifelse(is.na(B73_E) & !is.na(Mo17_E), "Mo17\ngenome\nonly", 
                                  ifelse(!is.na(B73_E) & is.na(Mo17_E), "B73\ngenome\nonly", 
                                         ifelse(B73_E != Mo17_E, "Discordant", "NA")))))) %>%
  mutate(A_Expression = ifelse(A == "Concordant", B73_A, 
                               ifelse(A == "B73\ngenome\nonly", B73_A, 
                                      ifelse(A == "Mo17\ngenome\nonly", Mo17_A, "NA")))) %>%
  mutate(B_Expression = ifelse(B == "Concordant", B73_B, 
                               ifelse(B == "B73\ngenome\nonly", B73_B, 
                                      ifelse(B == "Mo17\ngenome\nonly", Mo17_B, "NA")))) %>%
  mutate(C_Expression = ifelse(C == "Concordant", B73_C, 
                               ifelse(C == "B73\ngenome\nonly", B73_C, 
                                      ifelse(C == "Mo17\ngenome\nonly", Mo17_C, "NA")))) %>%
  mutate(D_Expression = ifelse(D == "Concordant", B73_D, 
                               ifelse(D == "B73\ngenome\nonly", B73_D, 
                                      ifelse(D == "Mo17\ngenome\nonly", Mo17_D, "NA")))) %>%
  mutate(E_Expression = ifelse(E == "Concordant", B73_E, 
                               ifelse(E == "B73\ngenome\nonly", B73_E, 
                                      ifelse(E == "Mo17\ngenome\nonly", Mo17_E, "NA"))))

deg_gene_types <- dplyr::select(deg_summary, Syntenic_Pair, B73_V4_GeneID, Mo17_CAU_GeneID, A:E) %>%
  gather("Section", "Gene_Type", 4:8)

deg_expression_patterns <- dplyr::select(deg_summary, Syntenic_Pair, B73_V4_GeneID, Mo17_CAU_GeneID, A_Expression:E_Expression) %>%
  gather("Section", "Expression_Type", 4:8) %>%
  mutate(Section = ifelse(Section == "A_Expression", "A", 
                          ifelse(Section == "B_Expression", "B", 
                                 ifelse(Section == "C_Expression", "C", 
                                        ifelse(Section == "D_Expression", "D", 
                                               ifelse(Section == "E_Expression", "E", "NA"))))))

deg_final_summary <- left_join(deg_gene_types, deg_expression_patterns, by = c("Syntenic_Pair", "B73_V4_GeneID", "Mo17_CAU_GeneID", "Section"))

### Determine the overlap of differentially expressed genes amongst the sections ###
#################################################################################################################
section_overlap <- filter(deg_final_summary, Gene_Type == "Concordant") %>% mutate(Gene_Type = 1) %>% spread(Section, Gene_Type)
section_overlap[is.na(section_overlap)] <- 0

### GO enrichment of genes differentially expressed between B73 & Mo17 at: all; husk-encased; and emerged section ###
#################################################################################################################
go_results_all <- GO_enrichment_func(myInterestingGenes = filter(section_overlap, A == "1", B == "1", C == "1", D == "1", E == "1")$B73_V4_GeneID, 
                                     go_file_path = "./data/V4_B73_maize_GAMER_GOs_Syntenic_Only.txt", 
                                     target_ontologies = "BP", p_value_cutoff = 0.001)

go_results_encased <- GO_enrichment_func(myInterestingGenes = filter(section_overlap, A == "1", B == "1", C == "1", D == "0", E == "0")$B73_V4_GeneID, 
                                         go_file_path = "./data/V4_B73_maize_GAMER_GOs_Syntenic_Only.txt", 
                                         target_ontologies = "BP", p_value_cutoff = 0.001)

go_results_emerged <- GO_enrichment_func(myInterestingGenes = filter(section_overlap, A == "0", B == "0", C == "0", D == "1", E == "1")$B73_V4_GeneID, 
                                         go_file_path = "./data/V4_B73_maize_GAMER_GOs_Syntenic_Only.txt", 
                                         target_ontologies = "BP", p_value_cutoff = 0.001)

### Graph results for each individual section ###
#################################################################################################################                    
### UpSet plot showing the overlap of differentially expressed genes across sections ###
upset(section_overlap, sets = c("A", "B", "C", "D", "E"),
      keep.order = TRUE, order.by = "freq", empty.intersections = "on",
      point.size = 6, line.size = 2.5, text.scale = 0, show.numbers = FALSE,
      mainbar.y.label = NULL)