### Load libraries ###
#################################################################################
library(tidyverse)
library(data.table)
library(DESeq2)
library(ggrepel)
library(grid)
library(scatterplot3d)

### Load data and compute stats ###
#################################################################################
library_info <- fread("./data/Library_Codes.csv")

syntenic_gene_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% 
  distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% 
  mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

b73_expressed_all <- fread("./data/RPKM_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, filter(library_info, Genotype == "B73")$Library_Name) %>%
  gather("Library_Name", "FPKM", 2:41) %>% 
  left_join(., dplyr::select(library_info, Library_Name, Section), by = "Library_Name") %>%
  group_by(GeneID, Section) %>% 
  summarise(FPKM = mean(FPKM, na.rm = TRUE)) %>% 
  filter(FPKM > 1)

mo17_expressed_all<- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, filter(library_info, Genotype == "Mo17")$Library_Name) %>%
  gather("Library_Name", "FPKM", 2:41) %>% 
  left_join(., dplyr::select(library_info, Library_Name, Section), by = "Library_Name") %>%
  group_by(GeneID, Section) %>% 
  summarise(FPKM = mean(FPKM, na.rm = TRUE)) %>% 
  filter(FPKM > 1)

b73_nonsyntenic_summary <- filter(b73_expressed_all, !GeneID %in% syntenic_gene_pairs$B73_V4_GeneID) %>% 
  dplyr::select(GeneID, Section) %>% 
  mutate(Genotype = "B73", Gene_Type = "Non-syntenic") %>% 
  group_by(Genotype, Section, Gene_Type) %>%
  summarise(Count = n())

mo17_nonsyntenic_summary <- filter(mo17_expressed_all, !GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID) %>% 
  dplyr::select(GeneID, Section) %>% 
  mutate(Genotype = "Mo17", Gene_Type = "Non-syntenic") %>% 
  group_by(Genotype, Section, Gene_Type) %>%
  summarise(Count = n())

b73_syntenic <- filter(b73_expressed_all, GeneID %in% syntenic_gene_pairs$B73_V4_GeneID) %>% 
  dplyr::select(GeneID, Section) %>%
  dplyr::rename("B73_V4_GeneID" = GeneID) %>% 
  inner_join(., syntenic_gene_pairs, by = "B73_V4_GeneID") %>% 
  ungroup() %>% 
  dplyr::select(Syntenic_Pair, Section) %>% 
  mutate(Expressed = "Yes", Genotype = "B73")

mo17_syntenic <- filter(mo17_expressed_all, GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID) %>% 
  dplyr::select(GeneID, Section) %>%
  dplyr::rename("Mo17_CAU_GeneID" = GeneID) %>% 
  inner_join(., syntenic_gene_pairs, by = "Mo17_CAU_GeneID") %>% 
  ungroup() %>% dplyr::select(Syntenic_Pair, Section) %>% 
  mutate(Expressed = "Yes", Genotype = "Mo17")

both_syntenic_summary <- bind_rows(b73_syntenic, mo17_syntenic) %>% spread(Genotype, Expressed) %>%
  mutate(Gene_Type = ifelse(is.na(B73), "Syntenic Mo17 Only", 
                            ifelse(is.na(Mo17), "Syntenic B73 Only", "Commonly Expressed Syntelog"))) %>%
  group_by(Section, Gene_Type) %>% summarise(Count = n())

expression_summary <- bind_rows(b73_nonsyntenic_summary, 
                                mo17_nonsyntenic_summary, 
                                filter(both_syntenic_summary, Gene_Type %in% c("Syntenic B73 Only", "Commonly Expressed Syntelog")) %>% mutate(Genotype = "B73"), 
                                filter(both_syntenic_summary, Gene_Type %in% c("Syntenic Mo17 Only", "Commonly Expressed Syntelog")) %>% mutate(Genotype = "Mo17")) %>%
  mutate(Gene_Type = ifelse(Gene_Type %in% c("Syntenic B73 Only", "Syntenic Mo17 Only"), "Genotype Specific Syntelog", Gene_Type)) %>%
  mutate(Gene_Type = factor(Gene_Type, levels = c("Non-syntenic", "Genotype Specific Syntelog", "Commonly Expressed Syntelog")))

b73_log2_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv")  %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name) 

mo17_log2_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name)

syntelog_log2_rpkm <- left_join(filter(b73_log2_rpkm, GeneID %in% syntenic_gene_pairs$B73_V4_GeneID), syntenic_gene_pairs, by = c("GeneID" = "B73_V4_GeneID")) %>%
  dplyr::select(-GeneID) %>%
  left_join(., (left_join(filter(mo17_log2_rpkm, GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID), syntenic_gene_pairs, by = c("GeneID" = "Mo17_CAU_GeneID")) %>% dplyr::select(-GeneID)), 
            by = "Syntenic_Pair") %>%
  dplyr::select(-B73_V4_GeneID, -Mo17_CAU_GeneID) %>% dplyr::select(Syntenic_Pair, everything()) %>% 
  drop_na() %>%
  mutate_at(vars(-Syntenic_Pair), .funs = funs(log2((. + 1))))
rownames(syntelog_log2_rpkm) <- NULL

# Principal component analysis
syntelog_log2_rpkm <- t(syntelog_log2_rpkm %>% column_to_rownames(var = "Syntenic_Pair"))
sample_order <- rownames(syntelog_log2_rpkm)

syntenic_pca <- prcomp(syntelog_log2_rpkm, center = TRUE)

syntenic_pca_results <- tibble(PC1 = syntenic_pca$x[, 1], 
                               PC2 = syntenic_pca$x[, 2], 
                               PC3 = syntenic_pca$x[, 3], 
                               Library_Name = sample_order) %>%
  left_join(., library_info, by = "Library_Name") %>%
  mutate(Genotype_Year = factor(paste(Genotype, ".", Year, sep = "")))

### Graph ###
#################################################################################
# Expression summary graph
p <- ggplot(expression_summary, aes(x = Genotype, y = Count, fill = Gene_Type)) + 
  facet_grid(. ~ Section) + 
  geom_bar(stat = "identity", color = "black", size = 1, width = 0.75) +
  theme_bw() +
  guides(fill = FALSE) +
  scale_fill_brewer(palette = "Greys") +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 36), 
        strip.background = element_rect(color = "black", size = 1.5),
        panel.grid.major.x = element_blank(), 
        strip.text = element_text(size = 52), 
        panel.border = element_rect(size = 3), 
        legend.position = "bottom") +
  labs(x = NULL, y = "# of Expressed Genes", fill = "Gene Type")

g <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl('strip-t', g$layout$name))

fills <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g)

### 3D PCA graphs ###
percentVars <- syntenic_pca$sdev^2/sum(syntenic_pca$sdev^2)
colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
colors <- colors[as.numeric(as.factor(syntenic_pca_results$Section))]
shapes = c(15, 22, 17, 24) 
shapes <- shapes[as.numeric(syntenic_pca_results$Genotype_Year)]

scatterplot3d(x = syntenic_pca_results$PC1, 
              y = syntenic_pca_results$PC3,
              z = syntenic_pca_results$PC2, 
              angle = 40, 
              color = colors, 
              pch = shapes, 
              cex.symbols = 4.5, 
              cex.axis = 2.75, 
              cex.lab = 2,
              label.tick.marks = TRUE,
              xlab = paste0("PC1: ", round(100 * percentVars[1], 0),"% variance"),
              ylab = paste0("PC3: ", round(100 * percentVars[3], 0),"% variance"), 
              zlab = paste0("PC2: ", round(100 * percentVars[2], 0),"% variance"))