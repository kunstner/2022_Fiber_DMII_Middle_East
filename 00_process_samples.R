
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)

library(ggpubr)
library(patchwork)

# Colors ------------------------------------------------------------------

# Colors for taxonomic abundance plots
taxNumberPlot <- 19
cols <- c(pals::tableau20(taxNumberPlot+1), "grey45")

colv <- c("steelblue", "orange")

# Data --------------------------------------------------------------------

ps_16s <- readRDS(file = "../P-003J_bac.phyloseq.RDS")
ps_its <- readRDS(file =  "../P-003J_fung.phyloseq.RDS")
metabo <- read_delim(file = "../Metabolites.csv", delim = ",")
cova   <- gdata::read.xls(xls = "../Metadata Metabelomics Characteristics.xlsx", sheet = 1)

metabo <- metabo %>% 
    dplyr::rename(Sample = "...1")
colnames(metabo) <- gsub(pattern = "\\.\\.\\.", replacement = ".", x = colnames(metabo))

cova <- cova %>% 
    dplyr::rename(Sample = Serial.No)

metabo$Sample

# Format phyloseq data ----------------------------------------------------

sample_names(ps_16s)
sample_names(ps_its)

# Reformat 16S rRNA data
# sample_data(ps_16s) %>% data.frame() %>% View()
tmp <- sample_data(ps_16s) %>% data.frame()
tmp <- tmp %>% 
    dplyr::rename(Sample = Your.Sample.ID)
tmp$Sample <- gsub(pattern = "C", replacement = "C0", x = tmp$Sample)
tmp$Sample <- gsub(pattern = "P", replacement = "P0", x = tmp$Sample)

tmp$Sample[tmp$Sample == "C03"] <- "C003"
tmp$Sample[tmp$Sample == "C04"] <- "C004"
tmp$Sample[tmp$Sample == "C05"] <- "C005"
tmp$Sample[tmp$Sample == "C06"] <- "C006"
tmp$Sample[tmp$Sample == "C07"] <- "C007"
tmp$Sample[tmp$Sample == "C08"] <- "C008"
tmp$Sample[tmp$Sample == "C09"] <- "C009"

tmp$Sample[tmp$Sample == "P03"] <- "P003"
tmp$Sample[tmp$Sample == "P04"] <- "P004"
tmp$Sample[tmp$Sample == "P05"] <- "P005"
tmp$Sample[tmp$Sample == "P06"] <- "P006"
tmp$Sample[tmp$Sample == "P07"] <- "P007"
tmp$Sample[tmp$Sample == "P08"] <- "P008"
tmp$Sample[tmp$Sample == "P09"] <- "P009"

rownames(tmp) <- tmp$Sample

sample_names(ps_16s) <- tmp$Sample
sample_data(ps_16s) <- tmp

# Reformat ITS data
# sample_data(ps_its) %>% data.frame() %>% View()
tmp <- sample_data(ps_its) %>% data.frame()
tmp <- tmp %>% 
    dplyr::rename(Sample = Your.Sample.ID)
tmp$Sample <- gsub(pattern = "C", replacement = "C0", x = tmp$Sample)
tmp$Sample <- gsub(pattern = "P", replacement = "P0", x = tmp$Sample)

tmp$Sample[tmp$Sample == "C03"] <- "C003"
tmp$Sample[tmp$Sample == "C04"] <- "C004"
tmp$Sample[tmp$Sample == "C05"] <- "C005"
tmp$Sample[tmp$Sample == "C06"] <- "C006"
tmp$Sample[tmp$Sample == "C07"] <- "C007"
tmp$Sample[tmp$Sample == "C08"] <- "C008"
tmp$Sample[tmp$Sample == "C09"] <- "C009"

tmp$Sample[tmp$Sample == "P03"] <- "P003"
tmp$Sample[tmp$Sample == "P04"] <- "P004"
tmp$Sample[tmp$Sample == "P05"] <- "P005"
tmp$Sample[tmp$Sample == "P06"] <- "P006"
tmp$Sample[tmp$Sample == "P07"] <- "P007"
tmp$Sample[tmp$Sample == "P08"] <- "P008"
tmp$Sample[tmp$Sample == "P09"] <- "P009"

rownames(tmp) <- tmp$Sample

sample_names(ps_its) <- tmp$Sample
sample_data(ps_its) <- tmp

# Select Samples ---------------------------------------------------------

tmp <- unique( c( sample_names(ps_16s), sample_names(ps_its), 
           cova$Sample, metabo$Sample ) )
tmp <- tmp[ tmp %in% sample_names(ps_16s) ]
tmp <- tmp[ tmp %in% cova$Sample ]
tmp <- tmp[ tmp %in% metabo$Sample ]

ps_16s <- subset_samples(physeq = ps_16s, Sample %in% tmp)
ps_its <- subset_samples(physeq = ps_its, Sample %in% tmp)
cova <- cova %>% dplyr::filter(Sample %in% tmp)
metabo <- metabo %>% dplyr::filter(Sample %in% tmp)

wilcox.test(Age ~ Status, data = cova)
boxplot(Age ~ Status, data = cova)

length(sample_names(ps_16s))
length(sample_names(ps_its))
nrow(cova)
nrow(metabo)

# add age
sample_names(ps_16s) == cova$Sample
sample_names(ps_its) == cova$Sample

sample_data(ps_16s)$Age <- cova$Age
sample_data(ps_its)$Age <- cova$Age

sample_data(ps_16s)$Weight <- cova$Weight.KG.
sample_data(ps_its)$Weight <- cova$Weight.KG.

sample_data(ps_16s)$BMI <- cova$BMI..Kg.m
sample_data(ps_its)$BMI <- cova$BMI..Kg.m

sample_data(ps_16s)$Gender <- cova$Gender
sample_data(ps_its)$Gender <- cova$Gender

# Screen -biome data ------------------------------------------------------

# 16S rRNA
tt <- tax_table(ps_16s) %>% data.frame()
table(tt$Phylum, useNA = 'always')
table(tt$Genus, useNA = 'always')
ps_16s <- subset_taxa(ps_16s, Phylum!="k__Bacteria_unclassified")

#ITS
tt <- tax_table(ps_its) %>% data.frame()
table(tt$Kingdom, useNA = 'always')
table(tt$Phylum, useNA = 'always')
ps_its <- subset_taxa(ps_its, Phylum!="Fungi_unclassified")

sample_data(ps_16s)$Shannon_16s <- estimate_richness(physeq = ps_16s, measures = "Shannon")$Shannon
sample_data(ps_its)$Shannon_ITS <- estimate_richness(physeq = ps_its, measures = "Shannon")$Shannon

sample_data(ps_16s)$Shannon_ITS <- sample_data(ps_its)$Shannon_ITS
sample_data(ps_its)$Shannon_16s <- sample_data(ps_16s)$Shannon_16s

# Save data for further analysis ------------------------------------------

saveRDS(object = ps_16s, file = "processed_ps_16s.RDS")
saveRDS(object = ps_its, file = "processed_ps_its.RDS")
saveRDS(object = metabo, file = "processed_metabo.RDS")
saveRDS(object = cova,   file = "processed_cova.RDS")

summary( sample_sums(ps_16s) )
summary( sample_sums(ps_its) )

