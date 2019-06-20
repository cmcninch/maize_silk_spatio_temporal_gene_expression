### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(readxl)
library(DESeq2)

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

### Compute Stats ###
#################################################################################################################
# Differential expression analysis between years
sections <- c("A", "B", "C", "D", "E")
degs_across_years <- ""
'%ni%' <- Negate('%in%')
for(i in 1:(length(sections))){
  b73_temp_df <- DiffExp_func(count_data = b73_counts, genotypes = "B73", sections = sections[i], design = ~ Year, contrast = c("Year", "2014", "2015")) %>% 
    mutate(Section = sections[i])
  mo17_temp_df <- DiffExp_func(count_data = mo17_counts, genotypes = "Mo17", sections = sections[i], design = ~ Year, contrast = c("Year", "2014", "2015")) %>% 
    mutate(Section = sections[i])
  
  b73_df <- left_join(dplyr::rename(b73_temp_df, "B73_V4_GeneID" = GeneID), syntenic_gene_pairs, by = "B73_V4_GeneID") %>% 
    mutate(Gene_Type = ifelse(Mo17_CAU_GeneID %in% mo17_temp_df$GeneID, "Commonly Differentially Expressed Syntelog", 
                              ifelse(Mo17_CAU_GeneID %ni% mo17_temp_df$GeneID & Mo17_CAU_GeneID %ni% syntenic_gene_pairs$Mo17_CAU_GeneID, "Non-syntenic", "Syntenic, genotype specific differential expression"))) %>%
    dplyr::select(Syntenic_Pair, B73_V4_GeneID, Mo17_CAU_GeneID, Section, Expression_Change, Gene_Type)
  
  mo17_df <- left_join(dplyr::rename(mo17_temp_df, "Mo17_CAU_GeneID" = GeneID), syntenic_gene_pairs, by = "Mo17_CAU_GeneID") %>% 
    mutate(Gene_Type = ifelse(B73_V4_GeneID %in% b73_temp_df$GeneID, "Commonly Differentially Expressed Syntelog", 
                              ifelse(B73_V4_GeneID %ni% b73_temp_df$GeneID & B73_V4_GeneID %ni% syntenic_gene_pairs$B73_V4_GeneID, "Non-syntenic", "Syntenic, genotype specific differential expression"))) %>%
    dplyr::select(Syntenic_Pair, B73_V4_GeneID, Mo17_CAU_GeneID, Section, Expression_Change, Gene_Type)
  
  syntenic_compare <- full_join(filter(b73_df, Gene_Type == "Commonly Differentially Expressed Syntelog") %>% dplyr::rename("B73_Expression_Change" = Expression_Change),
                                filter(mo17_df, Gene_Type == "Commonly Differentially Expressed Syntelog") %>% dplyr::rename("Mo17_Expression_Change" = Expression_Change), 
                                by = c("Syntenic_Pair", "B73_V4_GeneID", "Mo17_CAU_GeneID", "Gene_Type", "Section")) %>%
    mutate(Gene_Type = ifelse(B73_Expression_Change == Mo17_Expression_Change, "Syntenic, concordant differential expression", "Syntenic, discordant differential expression"))
  
  b73_df <- bind_rows((dplyr::rename(b73_df, "GeneID" = B73_V4_GeneID) %>% mutate(Genotype = "B73") %>% dplyr::select(GeneID, Genotype, Section, Expression_Change, Gene_Type) %>% filter(Gene_Type != "Commonly Differentially Expressed Syntelog")), 
                      (mutate(syntenic_compare, Genotype = "B73") %>% dplyr::rename("GeneID" = B73_V4_GeneID, "Expression_Change" = B73_Expression_Change) %>% dplyr::select(GeneID, Genotype, Section, Expression_Change, Gene_Type)))
  
  mo17_df <- bind_rows((dplyr::rename(mo17_df, "GeneID" = Mo17_CAU_GeneID) %>% mutate(Genotype = "Mo17") %>% dplyr::select(GeneID, Genotype, Section, Expression_Change, Gene_Type) %>% filter(Gene_Type != "Commonly Differentially Expressed Syntelog")), 
                       (mutate(syntenic_compare, Genotype = "Mo17") %>% dplyr::rename("GeneID" = Mo17_CAU_GeneID, "Expression_Change" = Mo17_Expression_Change) %>% dplyr::select(GeneID, Genotype, Section, Expression_Change, Gene_Type)))
  
  degs_across_years <- rbind(degs_across_years, b73_df, mo17_df)
}
degs_across_years_summary <- degs_across_years[-1, ] %>% group_by(Genotype, Section, Expression_Change, Gene_Type) %>% summarise(Gene_Count = n()) %>% 
  ungroup() %>% tidyr::complete(Genotype, Section, Gene_Type, Expression_Change, fill = list(Gene_Count = 0)) %>%
  mutate(Gene_Count = ifelse(Expression_Change == "Expression Increases", Gene_Count, -1*Gene_Count)) 

# Differential expression analysis across silk section transitions in individual years
sections <- c("A", "B", "C", "D", "E")
years <- c("2014", "2015")
degs_across_sections <- ""
'%ni%' <- Negate('%in%')
for(i in 1:length(years)){
  for(j in 1:(length(sections) - 1)){
    b73_temp_df <- DiffExp_func(count_data = b73_counts, genotypes = "B73", years = years[i], design = ~ Section, contrast = c("Section", sections[j], sections[j + 1])) %>% 
      mutate(Transition = paste(sections[j], "to", sections[j + 1], sep = " "))
    mo17_temp_df <- DiffExp_func(count_data = mo17_counts, genotypes = "Mo17", years = years[i], design = ~ Section, contrast = c("Section", sections[j], sections[j + 1])) %>% 
      mutate(Transition = paste(sections[j], "to", sections[j + 1], sep = " "))
    
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
                        (mutate(syntenic_compare, Genotype = "B73") %>% dplyr::rename("GeneID" = B73_V4_GeneID, "Expression_Change" = B73_Expression_Change) %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type))) %>% 
      mutate(Year = years[i])
    
    mo17_df <- bind_rows((dplyr::rename(mo17_df, "GeneID" = Mo17_CAU_GeneID) %>% mutate(Genotype = "Mo17") %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type) %>% filter(Gene_Type != "Commonly Differentially Expressed Syntelog")), 
                         (mutate(syntenic_compare, Genotype = "Mo17") %>% dplyr::rename("GeneID" = Mo17_CAU_GeneID, "Expression_Change" = Mo17_Expression_Change) %>% dplyr::select(GeneID, Genotype, Transition, Expression_Change, Gene_Type))) %>% 
      mutate(Year = years[i])
    
    degs_across_sections <- rbind(degs_across_sections, b73_df, mo17_df)
  }
}

degs_across_sections_summary <- degs_across_sections[-1, ] %>% group_by(Genotype, Year, Transition, Expression_Change, Gene_Type) %>% 
  summarise(Gene_Count = n()) %>% ungroup() %>%
  mutate(Gene_Count = ifelse(Expression_Change == "Expression Increases", Gene_Count, -1*Gene_Count))

### Graph ###
#################################################################################################################
# Differential expression across years in each section
ggplot(filter(degs_across_years_summary, Gene_Type != "Syntenic, discordant differential expression"), 
       aes(x = Section, y = Gene_Count, fill = reorder(Gene_Type, Gene_Count))) +
  guides(fill = FALSE) +
  geom_bar(stat = "identity", color = "black", size = 0.75, width = 0.75) +
  facet_grid(. ~ Genotype) +
  scale_fill_brewer(palette = "Greys") +
  geom_hline(yintercept = 0, size = 1.25, color = "black") +
  scale_y_continuous(limits = c(-1400, 800), breaks = c(-1250, -1000, -750, -500, -250,0, 250, 500, 750), 
                     labels = c("1250", "1000", "750", "500", "250", "0", "250", "500", "750")) +
  theme_bw() +
  theme(plot.title = element_blank(), 
        axis.title = element_text(size = 40), 
        axis.text = element_text(size = 32), 
        panel.border = element_rect(size = 5),
        legend.text = element_text(size = 16), 
        legend.background = element_rect(color = "black", size = 1), 
        legend.position = "right", 
        strip.background = element_rect(fill = "white", size = 2, color = "black"), 
        strip.text = element_text(size = 54, face = "bold"), 
        panel.spacing = unit(0.05, "lines")) +
  labs(y = "# of Differentially\nExpressed Genes", fill = NULL)

# Differential expression across section transitions
ggplot(filter(degs_across_sections_summary, Gene_Type != "Syntenic, discordant differential expression"), 
       aes(x = Transition, y = Gene_Count, fill = reorder(Gene_Type, Gene_Count))) +
  guides(fill = FALSE) +
  geom_bar(stat = "identity", color = "black", size = 0.75, width = 0.75) +
  facet_grid(Genotype ~ Year, scales = "free_x") +
  scale_fill_brewer(palette = "Greys") +
  geom_hline(yintercept = 0, size = 1.25, color = "black") +
  theme_bw() +
  theme(plot.title = element_blank(), 
        axis.title = element_text(size = 40), 
        axis.text = element_text(size = 32), 
        panel.border = element_rect(size = 5),
        legend.text = element_text(size = 16), 
        legend.background = element_rect(color = "black", size = 1), 
        legend.position = "right", 
        strip.background = element_rect(fill = "white", size = 2, color = "black"), 
        strip.text = element_text(size = 64, face = "bold"), 
        panel.spacing = unit(0.05, "lines")) +
  labs(y = "# of Differentially\nExpressed Genes", fill = NULL)