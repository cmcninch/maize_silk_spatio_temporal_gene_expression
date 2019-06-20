### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(DESeq2)
library(ComplexHeatmap)

### Load data ###
#################################################################################################################
library_codes <- fread("./data/Library_Codes.csv") %>%
                    mutate_at(vars(-Library_Name), factor)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
                    dplyr::select(GeneID, (filter(library_codes, Genotype == "B73"))$Library_Name)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
                      dplyr::select(GeneID, (filter(library_codes, Genotype == "Mo17"))$Library_Name)

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
                          distinct(B73_V4_GeneID, .keep_all = TRUE) %>% distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
                          arrange(B73_V4_GeneID) %>% mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

### Analyze data ###
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

sections <- c("A", "B", "C", "D", "E")
b73_degs_across_sections <- ""
for(i in 1:(length(sections) - 1)){
          temp_df <- DiffExp_func(count_data = b73_counts, genotypes = "B73", design = ~ Section + Year, contrast = c("Section", sections[i], sections[i + 1])) %>% 
                        mutate(Transition = paste(sections[i], "to", sections[i + 1], sep = " "))
          b73_degs_across_sections <- rbind(b73_degs_across_sections, temp_df)
}
b73_degs_across_sections <- b73_degs_across_sections[-1, ] %>% mutate(Genotype = "B73")

mo17_degs_across_sections <- ""
for(i in 1:(length(sections) - 1)){
          temp_df <- DiffExp_func(count_data = mo17_counts, genotypes = "Mo17", design = ~ Section + Year, contrast = c("Section", sections[i], sections[i + 1])) %>% 
                        mutate(Transition = paste(sections[i], "to", sections[i + 1], sep = " "))
          mo17_degs_across_sections <- rbind(mo17_degs_across_sections, temp_df)
}
mo17_degs_across_sections <- mo17_degs_across_sections[-1, ] %>% mutate(Genotype = "Mo17")

# Trends in expression across transitions all genes
b73_transition_trends <- dplyr::select(b73_degs_across_sections, GeneID, Transition, Expression_Change) %>%
                            spread(Transition, Expression_Change)
b73_transition_trends <- bind_rows(b73_transition_trends, 
                                   data.frame(GeneID = filter(b73_counts, !GeneID %in% b73_degs_across_sections$GeneID)$GeneID, 
                                              `A to B` = NA, `B to C` = NA, `C to D` = NA,
                                              `D to E` = NA, stringsAsFactors = FALSE, check.names = FALSE)) %>%
                            mutate(`A to B` = ifelse(is.na(`A to B`), 0, 
                                                     ifelse(`A to B` == "Expression Increases", 1, -1)), 
                                   `B to C` = ifelse(is.na(`B to C`), 0, 
                                                     ifelse(`B to C` == "Expression Increases", 1, -1)), 
                                   `C to D` = ifelse(is.na(`C to D`), 0, 
                                                     ifelse(`C to D` == "Expression Increases", 1, -1)), 
                                   `D to E` = ifelse(is.na(`D to E`), 0, 
                                                     ifelse(`D to E` == "Expression Increases", 1, -1))) %>%
                            group_by(`A to B`, `B to C`, `C to D`, `D to E`) %>% summarise(B73_Count = n())

mo17_transition_trends <- dplyr::select(mo17_degs_across_sections, GeneID, Transition, Expression_Change) %>%
                              spread(Transition, Expression_Change)
mo17_transition_trends <- bind_rows(mo17_transition_trends, 
                                    data.frame(GeneID = filter(mo17_counts, !GeneID %in% mo17_degs_across_sections$GeneID)$GeneID, 
                                               `A to B` = NA, `B to C` = NA, `C to D` = NA,
                                               `D to E` = NA, stringsAsFactors = FALSE, check.names = FALSE)) %>%
                              mutate(`A to B` = ifelse(is.na(`A to B`), 0, 
                                                       ifelse(`A to B` == "Expression Increases", 1, -1)), 
                                     `B to C` = ifelse(is.na(`B to C`), 0, 
                                                       ifelse(`B to C` == "Expression Increases", 1, -1)),
                                     `C to D` = ifelse(is.na(`C to D`), 0, 
                                                       ifelse(`C to D` == "Expression Increases", 1, -1)), 
                                     `D to E` = ifelse(is.na(`D to E`), 0, 
                                                       ifelse(`D to E` == "Expression Increases", 1, -1))) %>%
                              group_by(`A to B`, `B to C`, `C to D`, `D to E`) %>% summarise(Mo17_Count = n())

combined_transition_trends <- full_join(b73_transition_trends, mo17_transition_trends, by = c("A to B", "B to C", "C to D", "D to E"))
combined_transition_trends[is.na(combined_transition_trends)] <- 0
combined_transition_trends <- arrange(combined_transition_trends, -B73_Count, -Mo17_Count) %>%
                                  mutate(B73_Percentage = round(B73_Count/length(unique(b73_counts$GeneID)), digits = 3), 
                                         Mo17_Percentage = round(Mo17_Count/length(unique(mo17_counts$GeneID)), digits = 3))

# Trends in expression across transitions syntenic genes
b73_syntenic_transition_trends <- dplyr::select(b73_degs_across_sections, GeneID, Transition, Expression_Change) %>%
                                      spread(Transition, Expression_Change) %>% filter(GeneID %in% syntenic_gene_pairs$B73_V4_GeneID)
b73_syntenic_transition_trends <- bind_rows(b73_syntenic_transition_trends, 
                                            data.frame(GeneID = filter(b73_counts, !GeneID %in% b73_degs_across_sections$GeneID & GeneID %in% syntenic_gene_pairs$B73_V4_GeneID)$GeneID, 
                                                       `A to B` = NA, `B to C` = NA, `C to D` = NA,
                                                       `D to E` = NA, stringsAsFactors = FALSE, check.names = FALSE)) %>%
                                          mutate(`A to B` = ifelse(is.na(`A to B`), 0, 
                                                                   ifelse(`A to B` == "Expression Increases", 1, -1)), 
                                                 `B to C` = ifelse(is.na(`B to C`), 0, 
                                                                   ifelse(`B to C` == "Expression Increases", 1, -1)), 
                                                 `C to D` = ifelse(is.na(`C to D`), 0, 
                                                                   ifelse(`C to D` == "Expression Increases", 1, -1)), 
                                                 `D to E` = ifelse(is.na(`D to E`), 0, 
                                                                   ifelse(`D to E` == "Expression Increases", 1, -1))) %>%
                                          group_by(`A to B`, `B to C`, `C to D`, `D to E`) %>% summarise(B73_Count = n())

mo17_syntenic_transition_trends <- dplyr::select(mo17_degs_across_sections, GeneID, Transition, Expression_Change) %>%
                                        spread(Transition, Expression_Change) %>% filter(GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID)
mo17_syntenic_transition_trends <- bind_rows(mo17_syntenic_transition_trends, 
                                             data.frame(GeneID = filter(mo17_counts, !GeneID %in% mo17_degs_across_sections$GeneID & GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID)$GeneID, 
                                                        `A to B` = NA, `B to C` = NA, `C to D` = NA,
                                                        `D to E` = NA, stringsAsFactors = FALSE, check.names = FALSE)) %>%
                                        mutate(`A to B` = ifelse(is.na(`A to B`), 0, 
                                                                 ifelse(`A to B` == "Expression Increases", 1, -1)), 
                                               `B to C` = ifelse(is.na(`B to C`), 0, 
                                                                 ifelse(`B to C` == "Expression Increases", 1, -1)),
                                               `C to D` = ifelse(is.na(`C to D`), 0, 
                                                                 ifelse(`C to D` == "Expression Increases", 1, -1)), 
                                               `D to E` = ifelse(is.na(`D to E`), 0, 
                                                                 ifelse(`D to E` == "Expression Increases", 1, -1))) %>%
                                        group_by(`A to B`, `B to C`, `C to D`, `D to E`) %>% summarise(Mo17_Count = n())

combined_syntenic_transition_trends <- full_join(b73_syntenic_transition_trends, mo17_syntenic_transition_trends, by = c("A to B", "B to C", "C to D", "D to E"))
combined_syntenic_transition_trends[is.na(combined_syntenic_transition_trends)] <- 0
combined_syntenic_transition_trends <- arrange(combined_syntenic_transition_trends, -B73_Count, -Mo17_Count) %>%
                                          mutate(B73_Percentage = round(B73_Count/length(unique(syntenic_gene_pairs$B73_V4_GeneID)), digits = 3), 
                                                 Mo17_Percentage = round(Mo17_Count/length(unique(syntenic_gene_pairs$Mo17_CAU_GeneID)), digits = 3))

### Graph data ###
#################################################################################################################
# Used excel to color data frames. Red = 1; Grey = 0; Blue = -1