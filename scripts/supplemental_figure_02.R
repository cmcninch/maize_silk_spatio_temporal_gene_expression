### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(DESeq2)
library(VennDiagram)

### Load data ###
#################################################################################################################
library_codes <- fread("./data/Library_Codes.csv") %>%
                      mutate_at(vars(-Library_Name), factor)

b73_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv")  %>%
                dplyr::select(GeneID, (filter(library_codes, Genotype == "B73"))$Library_Name)

mo17_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
                dplyr::select(GeneID, (filter(library_codes, Genotype == "Mo17"))$Library_Name)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_codes, Genotype == "B73"))$Library_Name)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_codes, Genotype == "Mo17"))$Library_Name)

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
                        distinct(B73_V4_GeneID, .keep_all = TRUE) %>% distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
                        arrange(B73_V4_GeneID) %>% mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

### Analyze data ###
#################################################################################################################
### Section specific expression comparison B73 vs Mo17 using syntelogs ###
b73_expressed_gene_summary <- gather(b73_rpkm, "Library_Name", "RPKM", 2:41) %>% left_join(., library_codes, by = "Library_Name") %>%
                                group_by(GeneID, Genotype, Section) %>% summarise(RPKM = mean(RPKM, na.rm = TRUE)) %>%
                                mutate(Expressed = ifelse(RPKM > 1, 1, 0)) %>% ungroup() %>% dplyr::select(-RPKM) %>% 
                                spread(Section, Expressed) %>% mutate(`Any Section` = ifelse( A == 0 & B == 0 & C == 0 & D == 0 & E == 0, 0, 1))

mo17_expressed_gene_summary <- gather(mo17_rpkm, "Library_Name", "RPKM", 2:41) %>% left_join(., library_codes, by = "Library_Name") %>%
                                group_by(GeneID, Genotype, Section) %>% summarise(RPKM = mean(RPKM, na.rm = TRUE)) %>%
                                mutate(Expressed = ifelse(RPKM > 1, 1, 0)) %>% ungroup() %>% dplyr::select(-RPKM) %>% 
                                spread(Section, Expressed) %>% mutate(`Any Section` = ifelse( A == 0 & B == 0 & C == 0 & D == 0 & E == 0, 0, 1))

b73_syntenic_expressed_gene_summary <- gather(b73_rpkm, "Library_Name", "RPKM", 2:41) %>% left_join(., library_codes, by = "Library_Name") %>%
                                              group_by(GeneID, Genotype, Section) %>% summarise(RPKM = mean(RPKM, na.rm = TRUE)) %>%
                                              mutate(Expressed = ifelse(RPKM > 1, 1, 0)) %>% ungroup() %>% dplyr::select(-RPKM) %>% 
                                              filter(GeneID %in% syntenic_gene_pairs$B73_V4_GeneID) %>%
                                              spread(Section, Expressed) %>% mutate(`Any Section` = ifelse( A == 0 & B == 0 & C == 0 & D == 0 & E == 0, 0, 1))

mo17_syntenic_expressed_gene_summary <- gather(mo17_rpkm, "Library_Name", "RPKM", 2:41) %>% left_join(., library_codes, by = "Library_Name") %>%
                                              group_by(GeneID, Genotype, Section) %>% summarise(RPKM = mean(RPKM, na.rm = TRUE)) %>%
                                              mutate(Expressed = ifelse(RPKM > 1, 1, 0)) %>% ungroup() %>% dplyr::select(-RPKM) %>% 
                                              filter(GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID) %>%
                                              spread(Section, Expressed) %>% mutate(`Any Section` = ifelse( A == 0 & B == 0 & C == 0 & D == 0 & E == 0, 0, 1))

syntenic_overlap_b73_mo17 <- filter(syntenic_gene_pairs, 
                                    B73_V4_GeneID %in% filter(b73_expressed_gene_summary, `Any Section` == 1)$GeneID &
                                      Mo17_CAU_GeneID %in% filter(mo17_expressed_gene_summary, `Any Section` == 1)$GeneID)

### Graph data ###
#################################################################################################################
### Section specific expression ###
# Five way Venn diagram
venn_diagram_func <- function(df){
cat.just <- list(A = c(0.5, 0.70), B = c(-1, -0.75), C = c(1.5, -0.25), D = c(0.25, 0.5), E = c(1, 0.25))
grid.newpage()
draw.quintuple.venn(area1 = nrow(filter(df, A == 1)),
                    area2 = nrow(filter(df, B == 1)), 
                    area3 = nrow(filter(df, C == 1)), 
                    area4 = nrow(filter(df, D == 1)), 
                    area5 = nrow(filter(df, E == 1)), 
                    n12 = nrow(filter(df, A == 1 & B == 1)), 
                    n13 = nrow(filter(df, A == 1 & C == 1)), 
                    n14 = nrow(filter(df, A == 1 & D == 1)),
                    n15 = nrow(filter(df, A == 1 & E == 1)),
                    n23 = nrow(filter(df, B == 1 & C == 1)), 
                    n24 = nrow(filter(df, B == 1 & D == 1)), 
                    n25 = nrow(filter(df, B == 1 & E == 1)), 
                    n34 = nrow(filter(df, C == 1 & D == 1)), 
                    n35 = nrow(filter(df, C == 1 & E == 1)),
                    n45 = nrow(filter(df, D == 1 & E == 1)), 
                    n123 = nrow(filter(df, A == 1 & B == 1 & C == 1)), 
                    n124 = nrow(filter(df, A == 1 & B == 1 & D == 1)), 
                    n125 = nrow(filter(df, A == 1 & B == 1 & E == 1)), 
                    n134 = nrow(filter(df, A == 1 & C == 1 & D == 1)),
                    n135 = nrow(filter(df, A == 1 & C == 1 & E == 1)), 
                    n145 = nrow(filter(df, A == 1 & D == 1 & E == 1)),
                    n234 = nrow(filter(df, B == 1 & C == 1 & D == 1)), 
                    n235 = nrow(filter(df, B == 1 & C == 1 & E == 1)), 
                    n245 = nrow(filter(df, B == 1 & D == 1 & E == 1)), 
                    n345 = nrow(filter(df, C == 1 & D == 1 & E == 1)), 
                    n1234 = nrow(filter(df, A == 1 & B == 1 & C == 1 & D == 1)), 
                    n1235 = nrow(filter(df, A == 1 & B == 1 & C == 1 & E == 1)),
                    n1245 = nrow(filter(df, A == 1 & B == 1 & D == 1 & E == 1)), 
                    n1345 = nrow(filter(df, A == 1 & C == 1 & D == 1 & E == 1)), 
                    n2345 = nrow(filter(df, B == 1 & C == 1 & D == 1 & E == 1)), 
                    n12345 = nrow(filter(df, A == 1 & B == 1 & C == 1 & D == 1 & E == 1)), 
                    category = c("A", "B","C","D","E"),
                    lty = rep(1, 5), 
                    fill = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"), 
                    cex = 3, 
                    cat.just = cat.just,
                    cat.cex = 5)
}

venn_diagram_func(df = b73_expressed_gene_summary)
venn_diagram_func(df = mo17_expressed_gene_summary)
venn_diagram_func(df = b73_syntenic_expressed_gene_summary)
venn_diagram_func(df = mo17_syntenic_expressed_gene_summary)