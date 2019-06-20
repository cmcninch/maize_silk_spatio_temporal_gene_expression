### Load libraries ###
#################################################################################
library(tidyverse)
library(data.table)

### Load data ###
#################################################################################
photosynthesis_genes <- read_csv("./data/Photosynthesis_Genes.csv")

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% 
  distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% 
  mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

library_info <- fread("./data/Library_Codes.csv")

b73_log2_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name) %>%
  mutate_at(vars(-GeneID), .funs = funs(log2((. + 1)))) %>%
  gather("Library_Name", "Log2_RPKM", 2:ncol(.)) %>%
  left_join(., library_info, by = "Library_Name") %>%
  filter(GeneID %in% photosynthesis_genes$B73_GeneID)

mo17_log2_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name) %>%
  mutate_at(vars(-GeneID), .funs = funs(log2((. + 1)))) %>%
  gather("Library_Name", "Log2_RPKM", 2:ncol(.)) %>%
  left_join(., library_info, by = "Library_Name") %>%
  filter(GeneID %in% photosynthesis_genes$Mo17_GeneID)

b73_log2_rpkm_avg <- group_by(b73_log2_rpkm, GeneID, Section) %>%
  summarise(SE = sd(Log2_RPKM)/sqrt(n()),
            Log2_RPKM = mean(Log2_RPKM, na.rm = TRUE)) %>%
  left_join(., dplyr::select(photosynthesis_genes, B73_GeneID, Category, Order, Locus), by = c("GeneID" = "B73_GeneID")) %>%
  mutate(Genotype = "B73")

mo17_log2_rpkm_avg <- group_by(mo17_log2_rpkm, GeneID, Section) %>%
  summarise(SE = sd(Log2_RPKM)/sqrt(n()),
            Log2_RPKM = mean(Log2_RPKM, na.rm = TRUE)) %>%
  left_join(., dplyr::select(photosynthesis_genes, Mo17_GeneID, Category, Order, Locus), by = c("GeneID" = "Mo17_GeneID")) %>%
  mutate(Genotype = "Mo17")

avg_rpkm <- bind_rows(b73_log2_rpkm_avg, mo17_log2_rpkm_avg)

### Define functions ###
#################################################################################################################
line_graph_func <- function(category, num_cols = 1, scale_factor = 1) {
  
  formatter <- function(...){
    function(x) format(round(x, 1), ...)
  }
  
  print(ggplot(filter(avg_rpkm, Category == category),
               aes(x = Section, y = Log2_RPKM, color = Genotype, group = Genotype)) +
          geom_point(size = 6) +
          geom_line(size = 3) +
          guides(color = FALSE) +
          scale_x_discrete(expand = c(0.025, 0.025)) +
          scale_y_continuous(labels = formatter(nsmall = 1)) +
          geom_errorbar(aes(ymin = Log2_RPKM - SE, ymax = Log2_RPKM + SE), width = 0.1, size = 1.15 * scale_factor, color = "black") +
          facet_wrap(reorder(Locus, Order) ~ ., scales = "free_y", ncol = num_cols) +
          theme_bw() +
          theme(axis.text = element_text(size = 36), 
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                #axis.title.y = element_text(size = 48, margin = margin(0, 10, 0, 0)),
                strip.text = element_text(size = 48, margin = margin(0, 0, 5, 5), face = "italic"), 
                strip.background = element_blank(), 
                panel.border = element_rect(size = 4)) +
          labs(y = "Log2(RPKM+1)"))
}

### Compute stats ###
#################################################################################
t_test_ready_df <- bind_rows(left_join(b73_log2_rpkm, photosynthesis_genes, by = c("GeneID" = "B73_GeneID")) %>%
                               dplyr::select(GeneID, Genotype, Section, Locus, Log2_RPKM), 
                             left_join(mo17_log2_rpkm, photosynthesis_genes, by = c("GeneID" = "Mo17_GeneID")) %>%
                               dplyr::select(GeneID, Genotype, Section, Locus, Log2_RPKM))

t_test_results <- tibble(Locus = as.character(), Section = as.character(), P_Value = as.numeric())
loci <- unique(t_test_ready_df$Locus)
sections <- unique(t_test_ready_df$Section)
for(i in seq_along(loci)){
  for(j in seq_along(sections)){
    temp_result <- t.test(filter(t_test_ready_df, Locus == loci[i], Section == sections[j], Genotype == "B73")$Log2_RPKM, 
                          filter(t_test_ready_df, Locus == loci[i], Section == sections[j], Genotype == "Mo17")$Log2_RPKM)$p.value
    t_test_results <- bind_rows(t_test_results, 
                                tibble(Locus = loci[i], Section = sections[j], P_Value = temp_result))
  }
}

t_test_results <- mutate(t_test_results, Significant = ifelse(P_Value < 0.05, "Yes", "No"))

### Graph stats ###
#################################################################################
line_graph_func(category = "Light harvesting complex", num_cols = 2)
line_graph_func(category = "Photosystem I", num_cols = 2)
line_graph_func(category = "Photosystem II")
line_graph_func(category = "Cytb6f")
line_graph_func(category = "Ferredoxin")
line_graph_func(category = "Plastoquinone")