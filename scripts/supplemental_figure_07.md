#Plasma membrane and secreted peptide protocol
---
####Colton McNinch (June, 11th 2019)
* The following protocol highlights the individual steps taken to examine the number of C vs D DEGs localized in the plasma membrane or extracellular matrix. In addition, the number and type of Cysteine Rich Peptides are examined.

* The folowing protein sequences were used:

[B73 protein sequences]
(ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/pep/Zea_mays.AGPv4.pep.all.fa.gz)

[Mo17 protein sequences]
(https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/Zm00014a.proteins.fa.gz)

---
##Step 1: Obtain longest transcript sequence from C vs D DEGs.
#### Part A. Retrieve and unzip fasta files

```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/pep/Zea_mays.AGPv4.pep.all.fa.gz
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/Zm00014a.proteins.fa.gz

gunzip Zea_mays.AGPv4.pep.all.fa.gz
gunzip Zm00014a.proteins.fa.gz

```
#### Part B. Remove "*" from Mo17 file
```
tail -n+2 Zm00014a.proteins.fa | sed 's/*//g' Zm00014a.proteins.fa > Zm00014a.proteins.fa.new && mv Zm00014a.proteins.fa.new Zm00014a.proteins.fa
```
#### Part C. Retrieve `.fasta` files that contain all identified C to D differentially expressed genes by running the following Rscript `get_c_to_d_DEGs_fasta.R`

```
### Load libraries ###
#################################################################################################################
library(Biostrings)
library(tidyverse)
library(data.table)
library(DESeq2)

### Load data ###
#################################################################################################################
library_info <- fread("./data/Library_Codes.csv") %>%
  mutate_at(vars(-Library_Name), factor)

b73_counts <- fread("./data/Gene_Count_Matrix_with_B73_V4_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "B73"))$Library_Name)

mo17_counts <- fread("./data/Gene_Count_Matrix_with_CAU_Mo17_Genome.csv") %>%
  dplyr::select(GeneID, (filter(library_info, Genotype == "Mo17"))$Library_Name)
  
### Define functions ###
#################################################################################################################
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

### Differential expression analysis between C vs D sections ###
#################################################################################################################
b73_c_to_d_DEGs <- DiffExp_func(count_data = b73_counts, genotypes = "B73", design = ~ Section + Year, contrast = c("Section", "C", "D"))
mo17_c_to_d_DEGs <- DiffExp_func(count_data = mo17_counts, genotypes = "Mo17", design = ~ Section + Year, contrast = c("Section", "C", "D"))

### Create b73 and mo17 C vs D DEG fasta ###
#################################################################################################################
  b73_original_fasta <- readAAStringSet("./Zea_mays.AGPv4.pep.all.fa", seek.first.rec = TRUE)
  b73_longest_transcripts <- data.frame(Seq_Name = names(b73_original_fasta)) %>% 
                            mutate(Transcript_Name = substring(Seq_Name, regexpr("transcript:.+", Seq_Name)), 
                                   Gene_Name = substring(Seq_Name, regexpr("gene:.+", Seq_Name))) %>%
                            mutate(Transcript_Name = gsub(".*:","", gsub( " .*$", "", Transcript_Name)),
                                   Gene_Name = gsub(".*:","", gsub( " .*$", "", Gene_Name))) %>%
                            mutate(Length = b73_original_fasta@ranges@width) %>%
                            group_by(Gene_Name) %>% 
                            dplyr::slice(which.max(Length)) %>%
                            ungroup() %>%
                            filter(grepl("Zm", Gene_Name))
  
  b73_longest_transcripts <- as.data.frame(b73_longest_transcripts)
  b73_new_fasta <- b73_original_fasta[b73_longest_transcripts$Seq_Name]
  b73_new_fasta_given_names <- names(b73_new_fasta)
  names(b73_new_fasta) <- b73_longest_transcripts[match(b73_new_fasta_given_names, b73_longest_transcripts$Seq_Name), "Gene_Name"]
  b73_new_fasta <- b73_new_fasta[b73_c_to_d_DEGs$GeneID]
  writeXStringSet(b73_new_fasta, "./b73_c_to_d_DEGs.fa")
  
  mo17_original_fasta <- readAAStringSet("./Zm00014a.proteins.fa", seek.first.rec = TRUE)
  mo17_longest_transcripts <- data.frame(Seq_Name = names(mo17_original_fasta)) %>% 
                            mutate(Transcript_Name = Seq_Name, 
                                   Gene_Name = gsub( "_.*$", "", Seq_Name)) %>%
                            mutate(Length = mo17_original_fasta@ranges@width) %>%
                            group_by(Gene_Name) %>% 
                            dplyr::slice(which.max(Length)) %>%
                            filter(grepl("Zm", Gene_Name))
  
  mo17_longest_transcripts <- as.data.frame(mo17_longest_transcripts)
  mo17_new_fasta <- mo17_original_fasta[mo17_longest_transcripts$Seq_Name]
  mo17_new_fasta_given_names <- names(mo17_new_fasta)
  names(mo17_new_fasta) <- mo17_longest_transcripts[match(mo17_new_fasta_given_names, mo17_longest_transcripts$Seq_Name), "Gene_Name"]
  mo17_new_fasta <- mo17_new_fasta[mo17_c_to_d_DEGs$GeneID]
  writeXStringSet(mo17_new_fasta, "./mo17_c_to_d_DEGs.fa")
```
```
Rscript get_c_to_d_DEGs_fasta.R
```
---
## Step 2: Run SignalP 5.0 on all identified C to D differentially expressed genes to identify those with a signal peptide sequence. This can be done [here] (http://www.cbs.dtu.dk/services/SignalP/). Use the `Eukarya`, `Short output (no figures)` options.

---
## Step 3: Retrieve `.fasta` files that contain all identified C to D differentially expressed genes with signal sequences, with those sequences cleaved, by running the following Rscript `get_signal_c_to_d_DEGs_fasta.R` 
```
### Load libraries ###
#################################################################################################################
library(tidyverse)
library(Biostrings)

### Load data and get new fasta files ###
#################################################################################################################
b73_signal_c_to_d_DEGs <- read.delim("./b73_c_to_d_DEGs_SignalP.txt", comment.char = "#", sep = "", header = FALSE, 
                                     col.names = c("GeneID", "Prediction", "SP_Likelihood", "Other", as.character(1:2), "Cleavage_Site_Position", as.character(3:5))) %>%
  dplyr::select(GeneID, Prediction, SP_Likelihood, Cleavage_Site_Position) %>%
  mutate(Cleavage_Site_Position = as.numeric(substring(Cleavage_Site_Position, 1, 2)) + 1) %>%
  filter(Prediction != "OTHER")

mo17_signal_c_to_d_DEGs <- read.delim("./mo17_c_to_d_DEGs_SignalP.txt", comment.char = "#", sep = "", header = FALSE, 
                                      col.names = c("GeneID", "Prediction", "SP_Likelihood", "Other", as.character(1:2), "Cleavage_Site_Position", as.character(3:5))) %>%
  dplyr::select(GeneID, Prediction, SP_Likelihood, Cleavage_Site_Position) %>%
  mutate(Cleavage_Site_Position = as.numeric(substring(Cleavage_Site_Position, 1, 2)) + 1) %>%
  filter(Prediction != "OTHER")

b73_c_to_d_DEGs_fasta <- readAAStringSet("./b73_c_to_d_DEGs.fa", seek.first.rec = TRUE)
b73_signal_c_to_d_DEGs_fasta <- b73_c_to_d_DEGs_fasta[b73_signal_c_to_d_DEGs$GeneID]
b73_signal_c_to_d_DEGs_fasta <- padAndClip(b73_signal_c_to_d_DEGs_fasta,
                                           IRanges(start = b73_signal_c_to_d_DEGs$Cleavage_Site_Position, 
                                                   end = b73_signal_c_to_d_DEGs_fasta@ranges@width), 
                                           Lpadding.letter = "+", Rpadding.letter = "-", remove.out.of.view.strings = TRUE)

mo17_c_to_d_DEGs_fasta <- readAAStringSet("./mo17_c_to_d_DEGs.fa", seek.first.rec = TRUE)
mo17_signal_c_to_d_DEGs_fasta <- mo17_c_to_d_DEGs_fasta[mo17_signal_c_to_d_DEGs$GeneID]
mo17_signal_c_to_d_DEGs_fasta <- padAndClip(mo17_signal_c_to_d_DEGs_fasta,
                                           IRanges(start = mo17_signal_c_to_d_DEGs$Cleavage_Site_Position, 
                                                   end = mo17_signal_c_to_d_DEGs_fasta@ranges@width), 
                                           Lpadding.letter = "+", Rpadding.letter = "-", remove.out.of.view.strings = TRUE)

writeXStringSet(b73_signal_c_to_d_DEGs_fasta, "b73_cleaved_signal_c_to_d_DEGs.fa")
writeXStringSet(mo17_signal_c_to_d_DEGs_fasta, "mo17_cleaved_signal_c_to_d_DEGs.fa")

```
```
Rscript get_signal_c_to_d_DEGs_fasta.R
```
---
## Step 4: Run TMHMM 2.0 on all differentially expressed genes that have a signal peptide sequence [here] (http://www.cbs.dtu.dk/services/TMHMM/). Use the `One line per protein` option. Copy the output and save as a `.csv` file with no header.

---
## Step 5: Create a fasta file for proteins that don't have a transmembrane signal by running the following Rscript `get_nontransmembrane_fasta.R`.
```
### Load libraries ###
#################################################################################################################
library(tidyverse)
library(Biostrings)
library(data.table)

### Load data and get new fasta files ###
#################################################################################################################
### B73 ###
b73_c_to_d_DEGs_with_signal_sequence <- read.delim("./b73_cleaved_signal_c_to_d_DEGs_TMHMM.csv", sep = ",", header = FALSE, 
                                            col.names = c("GeneID", "Length", "Stat1", "Stat2", "Predicted_Helix", "Topology")) %>%
  dplyr::select(GeneID, Length, Predicted_Helix) %>%
  mutate(Length = as.numeric(gsub("len=", "", Length)), 
         Predicted_Helix = ifelse(Predicted_Helix == "PredHel=0", "No", "Yes"))

b73_non_transmembrane_secreted_c_to_d_DEGs <- filter(b73_c_to_d_DEGs_with_signal_sequence, Predicted_Helix == "No")

b73_cleaved_signal_c_to_d_DEGs_fa <- readAAStringSet("b73_cleaved_signal_c_to_d_DEGs.fa", seek.first.rec = TRUE)

b73_non_transmembrane_secreted_c_to_d_DEGs_fa <- b73_cleaved_signal_c_to_d_DEGs_fa[b73_non_transmembrane_secreted_c_to_d_DEGs$GeneID]

writeXStringSet(b73_non_transmembrane_secreted_c_to_d_DEGs_fa, "./b73_nontransmembrane_signal_c_to_d_DEGs.fa")

### Mo17 ###
mo17_c_to_d_DEGs_with_signal_sequence <- read.delim("./mo17_cleaved_signal_c_to_d_DEGs_TMHMM.csv", sep = ",", header = FALSE, 
                                                    col.names = c("GeneID", "Length", "Stat1", "Stat2", "Predicted_Helix", "Topology")) %>%
  dplyr::select(GeneID, Length, Predicted_Helix) %>%
  mutate(Length = as.numeric(gsub("len=", "", Length)), 
         Predicted_Helix = ifelse(Predicted_Helix == "PredHel=0", "No", "Yes"))
         
mo17_non_transmembrane_secreted_c_to_d_DEGs <- filter(mo17_c_to_d_DEGs_with_signal_sequence, Predicted_Helix == "No")

mo17_cleaved_signal_c_to_d_DEGs_fa <- readAAStringSet("mo17_cleaved_signal_c_to_d_DEGs.fa", seek.first.rec = TRUE)

mo17_non_transmembrane_secreted_c_to_d_DEGs_fa <- mo17_cleaved_signal_c_to_d_DEGs_fa[mo17_non_transmembrane_secreted_c_to_d_DEGs$GeneID]

writeXStringSet(mo17_non_transmembrane_secreted_c_to_d_DEGs_fa, "./mo17_nontransmembrane_signal_c_to_d_DEGs.fa")
```
---
## Step 6: Scan non-transmembrane proteins that had a signal sequence for Cysteine Rich Peptide (CRP) motifs. This will be done with hmmer 3.2 scans using previously made [CRP Hidden Markov Models (HMMs)](https://doi.org/10.1111/j.1365-313X.2007.03136.x) 
```
hmmsearch --tblout b73_CRPs.tab -E 1e-5 plant_CRPs.hmm b73_nontransmembrane_signal_c_to_d_DEGs.fa
hmmsearch --tblout mo17_CRPs.tab -E 1e-5 plant_CRPs.hmm mo17_nontransmembrane_signal_c_to_d_DEGs.fa
```
---
## Step 7: Summarise results and make figures by running the Rscript `summarise_results.r`
```
### Load libraries ###
#################################################################################################################
library(tidyverse)
library(Biostrings)
library(readxl)
library(data.table)
library(ggforce)

### Load data and format for consolidated summaries ###
#################################################################################################################
### B73 total summaries ###
b73_signal_scans <- read.delim("./b73_c_to_d_DEGs_SignalP.txt", comment.char = "#", sep = "", header = FALSE, 
                             col.names = c("GeneID", "Prediction", "SP_Likelihood", "Other", as.character(1:2), "Cleavage_Site_Position", as.character(3:5)))

b73_non_signal <- filter(b73_signal_scans, Prediction == "OTHER") %>% 
  dplyr::select(GeneID) %>%
  mutate(Gene_Type = "No Signal Sequence")

b73_transmembrane <- read.delim("./b73_cleaved_signal_c_to_d_DEGs_TMHMM.csv", sep = ",", header = FALSE, 
                                col.names = c("GeneID", "Length", "Stat1", "Stat2", "Predicted_Helix", "Topology")) %>%
  filter(Predicted_Helix != "PredHel=0") %>%
  dplyr::select(GeneID) %>%
  mutate(Gene_Type = "Transmembrane")

b73_crps <- read.delim("./b73_CRPs.tab", comment.char = "#", sep = "", header = FALSE, 
                       col.names = c("GeneID", "1", "CRP_Class", "2", "E_Value", "3", as.character(4:16))) %>%
  dplyr::select(GeneID, E_Value, CRP_Class) %>%
  group_by(GeneID) %>%
  dplyr::slice(which.min(E_Value)) %>%
  dplyr::select(GeneID) %>%
  mutate(Gene_Type = "CRP")

b73_non_crps_peps <- tibble(GeneID = setdiff(names(readAAStringSet("./b73_nontransmembrane_signal_c_to_d_DEGs.fa", seek.first.rec = TRUE)), b73_crps$GeneID)) %>%
  mutate(Gene_Type = "Non-Cys-Rich Peptide")

b73_c_to_d_DEGs <- bind_rows(b73_non_signal, b73_transmembrane, b73_crps, b73_non_crps_peps) %>%
  mutate(Genotype = "B73")

### Mo17 total summaries ###
mo17_signal_scans <- read.delim("./mo17_c_to_d_DEGs_SignalP.txt", comment.char = "#", sep = "", header = FALSE, 
                               col.names = c("GeneID", "Prediction", "SP_Likelihood", "Other", as.character(1:2), "Cleavage_Site_Position", as.character(3:5)))

mo17_non_signal <- filter(mo17_signal_scans, Prediction == "OTHER") %>% 
  dplyr::select(GeneID) %>%
  mutate(Gene_Type = "No Signal Sequence")

mo17_transmembrane <- read.delim("./mo17_cleaved_signal_c_to_d_DEGs_TMHMM.csv", sep = ",", header = FALSE, 
                                col.names = c("GeneID", "Length", "Stat1", "Stat2", "Predicted_Helix", "Topology")) %>%
  filter(Predicted_Helix != "PredHel=0") %>%
  dplyr::select(GeneID) %>%
  mutate(Gene_Type = "Transmembrane")

mo17_crps <- read.delim("./mo17_CRPs.tab", comment.char = "#", sep = "", header = FALSE, 
                       col.names = c("GeneID", "1", "CRP_Class", "2", "E_Value", "3", as.character(4:16))) %>%
  dplyr::select(GeneID, E_Value, CRP_Class) %>%
  group_by(GeneID) %>%
  dplyr::slice(which.min(E_Value)) %>%
  dplyr::select(GeneID) %>%
  mutate(Gene_Type = "CRP")

mo17_non_crps_peps <- tibble(GeneID = setdiff(names(readAAStringSet("./mo17_nontransmembrane_signal_c_to_d_DEGs.fa", seek.first.rec = TRUE)), mo17_crps$GeneID)) %>%
  mutate(Gene_Type = "Non-Cys-Rich Peptide")

mo17_c_to_d_DEGs <- bind_rows(mo17_non_signal, mo17_transmembrane, mo17_crps, mo17_non_crps_peps) %>%
  mutate(Genotype = "Mo17")

c_to_d_DEGs <- bind_rows(b73_c_to_d_DEGs, mo17_c_to_d_DEGs) %>% 
  group_by(Gene_Type, Genotype) %>%
  summarise(Count = n())

### Load data and format for CRP summaries ###
#################################################################################################################
crp_info <- read_xls("./tpj_3136_sm_suppmat/TPJ3136TableS1.xls", skip = 1)
crp_info <- crp_info[1:2]
colnames(crp_info) <- c("CRP_ClassID", "CRP_Class")

### B73 CRP summaries ###
b73_crps <- read.delim("./b73_CRPs.tab", comment.char = "#", sep = "", header = FALSE, 
                       col.names = c("GeneID", "1", "CRP_ClassID", "2", "E_Value", "3", as.character(4:16))) %>%
  dplyr::select(GeneID, E_Value, CRP_ClassID) %>%
  group_by(GeneID) %>%
  dplyr::slice(which.min(E_Value)) %>%
  dplyr::select(GeneID, CRP_ClassID) %>%
  mutate(CRP_ClassID = gsub(".trim", "", CRP_ClassID), 
         Genotype = "B73")

### Mo17 CRP summaries ###
mo17_crps <- read.delim("./mo17_CRPs.tab", comment.char = "#", sep = "", header = FALSE, 
                       col.names = c("GeneID", "1", "CRP_ClassID", "2", "E_Value", "3", as.character(4:16))) %>%
  dplyr::select(GeneID, E_Value, CRP_ClassID) %>%
  group_by(GeneID) %>%
  dplyr::slice(which.min(E_Value)) %>%
  dplyr::select(GeneID, CRP_ClassID) %>%
  mutate(CRP_ClassID = gsub(".trim", "", CRP_ClassID), 
         Genotype = "Mo17")

### Combined CRP summaries ###
crps <- bind_rows(b73_crps, mo17_crps) %>% left_join(., crp_info, by = "CRP_ClassID") %>%
  mutate(CRP_Class = ifelse(CRP_Class %in% c("Chitinase (Class I/II) / Hevein",
                                             "Chitinase (Class IV) / Hevein", 
                                             "Chitinase / Hevein / PR-4 / Wheatwin2"), "Chitinase / Hevein",
                            ifelse(CRP_Class == "Lipid transfer protein (BETL4)", "Lipid transfer protein", 
                            ifelse(CRP_Class %in% c("Kunitz type trypsin inhibitor", 
                                                    "Novel GRP family", 
                                                    "STIG1", 
                                                    "Proteinase inhibitor II", 
                                                    "GASA/GAST/Snakin", 
                                                    "LCR", 
                                                    "Thionin related"), "Other", 
                            ifelse(CRP_Class == "pollen Ole e I family allergens", "Pollen allergen Ole e I family", CRP_Class))))) %>%
  group_by(Genotype, CRP_Class) %>%
  summarise(Count = n())

### Overlap of syntenic CRPs ###
syntenic_gene_pairs <- fread("~/Box Sync/McNinch Projects/Slice&Dice_RNA_Seq_Manuscript/data/V4_B73_vs_CAU_Mo17_Syntenic_Orthologs.csv") %>%
  distinct(B73_V4_GeneID, .keep_all = TRUE) %>% 
  distinct(Mo17_CAU_GeneID, .keep_all = TRUE) %>%
  arrange(B73_V4_GeneID) %>% 
  mutate(Syntenic_Pair = paste("Gene Pair", seq(1:nrow(.))))

non_syntenic_crps <- bind_rows(filter(b73_crps, !GeneID %in% syntenic_gene_pairs$B73_V4_GeneID), 
                               filter(mo17_crps, !GeneID %in% syntenic_gene_pairs$Mo17_CAU_GeneID)) %>%
  mutate(Gene_Type = "Non-syntenic") %>%
  dplyr::select(GeneID, Genotype, Gene_Type)

b73_syntenic_crps <- left_join(b73_crps, syntenic_gene_pairs, by = c("GeneID" = "B73_V4_GeneID")) %>%
  drop_na()

mo17_syntenic_crps <- left_join(mo17_crps, syntenic_gene_pairs, by = c("GeneID" = "Mo17_CAU_GeneID")) %>%
  drop_na()
  
b73_syntenic_crps <- mutate(b73_syntenic_crps, Gene_Type = ifelse(Mo17_CAU_GeneID %in% mo17_syntenic_crps$GeneID, 
                                                      "Commonly Differentially Expressed Syntelog", "Syntenic, genotype specific differential expression"))
  
mo17_syntenic_crps <- mutate(mo17_syntenic_crps, Gene_Type = ifelse(B73_V4_GeneID %in% b73_syntenic_crps$GeneID, 
                                                      "Commonly Differentially Expressed Syntelog", "Syntenic, genotype specific differential expression"))

syntenic_crps <- bind_rows(b73_syntenic_crps, mo17_syntenic_crps) %>%
  dplyr::select(GeneID, Genotype, Gene_Type)

crp_overlap <- bind_rows(non_syntenic_crps, syntenic_crps) %>% group_by(Genotype, Gene_Type) %>%
  summarise(Count = n()) %>%
  mutate(Gene_Type = factor(Gene_Type, levels = c("Non-syntenic", 
                                                  "Syntenic, genotype specific differential expression", 
                                                  "Commonly Differentially Expressed Syntelog")))

### Create figures ###
#################################################################################################################
### Overall summaries ###
ggplot(data = c_to_d_DEGs, aes(x = factor(1), y = Count, fill = Gene_Type)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") +
  coord_polar(theta = "y") +
  facet_grid(. ~ Genotype) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(), 
        panel.border = element_blank(),
        legend.text = element_text(size = 20), 
        legend.key = element_rect(color = "black"), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 40)) +
  labs(x = NULL, y = NULL, fill = NULL)

### CRP summaries ###
ggplot(data = crps, aes(x = factor(1), y = Count, fill = reorder(CRP_Class, Count))) + 
  geom_bar(stat = "identity", position = "fill", color = "black") +
  coord_polar(theta = "y") +
  facet_grid(. ~ Genotype) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(), 
        panel.border = element_blank(),
        legend.text = element_text(size = 20), 
        legend.key = element_rect(color = "black"), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 40)) +
  labs(x = NULL, y = NULL, fill = NULL)

### CRP overlap ###
ggplot(crp_overlap, aes(x = Genotype, y = Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.75, width = 0.75) +
  scale_fill_brewer(palette = "Greys") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(plot.title = element_blank(), 
        axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), 
        panel.border = element_rect(size = 5),
        legend.text = element_text(size = 16), 
        legend.background = element_rect(color = "black", size = 1)) +
  labs(y = "# of Differentially\nExpressed CRPs", fill = NULL, x = NULL)

```
```
Rscript summarise_results.R
```