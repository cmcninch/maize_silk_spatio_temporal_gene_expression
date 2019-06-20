### Load libraries ###
#################################################################################################################
library(tidyverse)
library(data.table)
library(stats)
library(ggpubr)
library(broom)

### Load data ###
#################################################################################################################
library_codes <- fread("./data/Library_Codes.csv") %>%
  dplyr::select(Library_Name, Genotype)

b73_v4_align_stats <- fread("./data/Alignment_Summary_with_B73_V4_Genome.csv") %>%
  mutate(Genome = "B73\nGenome") %>% dplyr::select(Library_Name, Genome, everything())

cau_mo17_align_stats <- fread("./data/Alignment_Summary_with_CAU_Mo17_Genome.csv") %>%
  mutate(Genome = "Mo17\nGenome") %>% dplyr::select(Library_Name, Genome, everything())

b73_rpkm <- fread("./data/RPKM_with_B73_V4_Genome.csv")  %>%
  dplyr::select(GeneID, (filter(library_codes, Genotype == "B73"))$Library_Name)

mo17_rpkm <- fread("./data/RPKM_with_CAU_Mo17_genome.csv") %>%
  dplyr::select(GeneID, (filter(library_codes, Genotype == "Mo17"))$Library_Name)

all_rpkm <- bind_rows(b73_rpkm, mo17_rpkm) %>% gather("Library_Name", "RPKM", 2:81)
all_rpkm <- left_join(all_rpkm, fread("./data/Library_Codes.csv"), by = c("Library_Name"))

syntenic_pairs <- fread("./data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv")

### Compute stats ###
#################################################################################################################
# Unique and Unaligned percentages ANOVA tests
combined_align_stats <- bind_rows(b73_v4_align_stats, cau_mo17_align_stats) %>%
  inner_join(., library_codes, by = "Library_Name") %>%
  mutate(Unique_Percent = Uniquely_Aligned_Reads/Total_Reads, 
         Unaligned_Percent = Unaligned_Reads/Total_Reads) %>%
  dplyr::select(Library_Name, Genome, Genotype, Unique_Percent, Unaligned_Percent) %>%
  gather("Stat", "Percent", 4:5) %>%
  mutate(Genotype = ifelse(Genotype == "B73", "B73 Samples", "Mo17 Samples"))


anova_compute_func <- function(stat, genotype){
  AOV <- aov(Percent ~ Genome, data = filter(combined_align_stats, Stat == stat, Genotype == genotype))
  tidy_aov <- tidy(AOV)
  
  p_value <- formatC(tidy_aov$p.value[1], format = "e", digits = 2)
  r_squared <- round(tidy_aov$sumsq[1]/(tidy_aov$sumsq[1] + tidy_aov$sumsq[2]), digits = 2)
  
  anova_results <- tibble(Genotype = genotype, Stat = stat, P_Value = p_value, R_Squared = r_squared)
  return(anova_results)
}

combined_align_anova_stats = bind_rows(anova_compute_func(stat = "Unaligned_Percent", genotype = "B73 Samples"), 
                                       anova_compute_func(stat = "Unaligned_Percent", genotype = "Mo17 Samples"),
                                       anova_compute_func(stat = "Unique_Percent", genotype = "B73 Samples"), 
                                       anova_compute_func(stat = "Unique_Percent", genotype = "Mo17 Samples"))

# Expression summary all genes
summary_df <- data.frame(Genotype = as.character(), Section = as.character(), SE = as.numeric(), 
                         Mean_Expressed_Genes = as.numeric(), Threshold = as.character(), stringsAsFactors = FALSE)
thresholds <- as.character(c(0, 1, 10, 100, 1000, 10000))
for(i in 1:(length(thresholds) - 1)){
  temp_df <- dplyr::select(all_rpkm, -Library_Name, -Location_Of_Silk) %>% group_by(Genotype, Section, Year, Rep) %>%
    summarise(Expressed_Genes = sum(RPKM >= as.numeric(thresholds[i]) & RPKM < as.numeric(thresholds[i+1]), na.rm = TRUE)) %>%
    group_by(Genotype, Section, Year) %>%
    summarise(SE = sd(Expressed_Genes)/sqrt(n()), 
              Mean_Expressed_Genes = mean(Expressed_Genes)) %>%
    mutate(Threshold = paste(thresholds[i], "-", thresholds[i+1], sep = ""))
  summary_df <- bind_rows(summary_df, temp_df)
}
summary_df <- summary_df %>% mutate(Threshold = factor(Threshold, levels = c("0-1", "1-10", "10-100", "100-1000", "1000-10000"))) %>%
  mutate(Genotype_Year = paste(Genotype, "-", Year, sep = ""))

# Expression summary syntenic genes
syntenic_summary_df <- data.frame(Genotype = as.character(), Section = as.character(), SE = as.numeric(), 
                                  Mean_Expressed_Genes = as.numeric(), Threshold = as.character(), stringsAsFactors = FALSE)
thresholds <- as.character(c(0, 1, 10, 100, 1000, 10000))
syntenic_rpkm <- filter(all_rpkm, GeneID %in% syntenic_pairs$B73_V4_GeneID | GeneID %in% syntenic_pairs$Mo17_CAU_GeneID)
for(i in 1:(length(thresholds) - 1)){
  temp_df <- dplyr::select(syntenic_rpkm, -Library_Name, -Location_Of_Silk) %>% group_by(Genotype, Section, Year, Rep) %>%
    summarise(Expressed_Genes = sum(RPKM >= as.numeric(thresholds[i]) & RPKM < as.numeric(thresholds[i+1]), na.rm = TRUE)) %>%
    group_by(Genotype, Section, Year) %>%
    summarise(SE = sd(Expressed_Genes)/sqrt(n()), 
              Mean_Expressed_Genes = mean(Expressed_Genes)) %>%
    mutate(Threshold = paste(thresholds[i], "-", thresholds[i+1], sep = ""))
  syntenic_summary_df <- bind_rows(syntenic_summary_df, temp_df)
}
syntenic_summary_df <- syntenic_summary_df %>% mutate(Threshold = factor(Threshold, levels = c("0-1", "1-10", "10-100", "100-1000", "1000-10000"))) %>%
  mutate(Genotype_Year = paste(Genotype, "-", Year, sep = ""))


### Create figure ###
#################################################################################################################
### Boxplots ###
plot_func <- function(stat, trim_on, y_lims, y_breaks, y_labs, y_adj = 1, x_lab_size, x_tick_size, y_axis_title, 
                      x_stat_label_pos = c(1.5, 1.5), y_stat_label_pos = c(0.9, 0.9)){
  ggplot(filter(combined_align_stats, Stat == stat), aes(x = Genome, y = Percent)) +
    geom_boxplot(width = 0.75, size = 1.25, color = "black", aes(fill = Genome)) +
    facet_grid(. ~ Genotype) + 
    theme_bw() + guides(fill = FALSE) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks, labels = y_labs) +
    theme(axis.text.x = element_text(size = x_lab_size), 
          axis.text.y = element_text(size = 32, margin = margin(t = 0, r = 0, b = 0, l = 20)),
          axis.ticks.x = element_line(size = x_tick_size),
          axis.text = element_text(size = 32),
          axis.title = element_text(size = 32),
          strip.text = element_text(size = 40), 
          strip.background = element_rect(fill = "orange", size = 2), 
          panel.border = element_rect(size = 4)) +
    labs(x = NULL, y = y_axis_title)
}

plot_func(stat = "Unique_Percent", y_lims = c(0.5, 0.95), y_breaks = seq(0.5, 0.9, by = 0.1),
          x_stat_label_pos = c(1.125, 1.125), y_stat_label_pos = c(0.9, 0.9),
          y_labs = c("50%", "60%", "70%", "80%", "90%"), x_lab_size = 40, x_tick_size = 0,
          y_axis_title = "% of Reads with Unique Alignments")

plot_func(stat = "Unaligned_Percent", y_lims = c(0.05, 0.415), y_breaks = seq(0.1, 0.4, by = 0.1), 
          x_stat_label_pos = c(1.125, 1.125), y_stat_label_pos = c(0.375, 0.375),
          y_labs = c("10%", "20%", "30%", "40%"), y_adj = 1.155, x_lab_size = 40, x_tick_size = 1,
          y_axis_title = "% of Reads with no Alignments")

### Dotplot ###
# Expression trends graph all genes
ggplot(summary_df, aes(x = Threshold, y = Mean_Expressed_Genes, color = Genotype_Year, group = Genotype_Year)) +
  facet_grid(. ~ Section) +
  geom_point(size = 5, position = position_dodge(width = 0.35)) +
  guides(color = guide_legend(override.aes = list(size = 14))) +
  scale_y_continuous(breaks = c(1000, 5000, 10000, 15000, 20000)) +
  theme_bw() +
  theme(axis.title = element_text(size = 32), 
        axis.text.x = element_text(size = 32, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 32),
        strip.text = element_text(size = 40),
        legend.text = element_text(size = 32),
        plot.title = element_text(hjust = 0.5, size = 48),
        panel.border = element_rect(size = 4), 
        strip.background = element_rect(size = 2.5, fill = "orange")) +
  labs(x = "RPKM", y = "# of Genes ", color = NULL, title = "Expression Trends of all Genes")

# Expression trends graph syntenic genes
ggplot(syntenic_summary_df, aes(x = Threshold, y = Mean_Expressed_Genes, color = Genotype_Year, group = Genotype_Year)) +
  facet_grid(. ~ Section) +
  geom_point(size = 5, position = position_dodge(width = 0.35)) +
  guides(color = guide_legend(override.aes = list(size = 14))) +
  scale_y_continuous(limits = c(0, 10000), breaks = c(2500, 5000, 7500, 10000)) +
  theme_bw() +
  theme(axis.title = element_text(size = 32), 
        axis.text.x = element_text(size = 32, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 32),
        strip.text = element_text(size = 40),
        legend.text = element_text(size = 32),
        plot.title = element_text(hjust = 0.5, size = 48),
        panel.border = element_rect(size = 4), 
        strip.background = element_rect(size = 2.5, fill = "orange")) +
  labs(x = "RPKM", y = "# of Genes ", color = NULL, title = "Expression Trends of Syntenic Genes")