# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)

library(phyloseq)

# library(caret)
library(corrplot)

# Colors ------------------------------------------------------------------

# Colors for taxonomic abundance plots
taxNumberPlot <- 19
cols <- c(pals::tableau20(taxNumberPlot+1), "grey85")

colv <- c("steelblue", "orange")

seed <- 42
set.seed(seed)

# Get data ------------------------------------------------------

ps_16s <- readRDS(file = "data/processed_ps_16s.RDS")
ps_its <- readRDS(file = "data/processed_ps_its.RDS")
metabo <- readRDS(file = "data/processed_metabo.RDS")
cova   <- readRDS(file = "data/processed_cova.RDS")

ps_16s <- ps_16s %>% tax_glom(physeq = ., taxrank = "Genus")
ps_its <- ps_its %>% tax_glom(physeq = ., taxrank = "Genus")

bact <- otu_table(object = ps_16s) %>% data.frame()
fung <- otu_table(object = ps_its) %>% t() %>% data.frame()

rownames(bact) == rownames(fung) 

metabo <- metabo %>% 
    column_to_rownames("Sample") 

# remove taxa with low variability
bact <- bact[ , apply(bact, 2, var) > 10 ]
fung <- fung[ , apply(fung, 2, var) > 10 ]
metabo <- metabo[ , apply(metabo, 2, var) > 10 ]

ps_16s <- prune_taxa(taxa = colnames(bact), x = ps_16s) 
ps_its <- prune_taxa(taxa = colnames(fung), x = ps_its) 

# Normalise with CLR

norm_b.clr <- microbiome::transform(x = ps_16s, transform = "clr") %>% 
    otu_table() %>% data.frame()
colnames(norm_b.clr) <- data.frame( tax_table( object = ps_16s ) )[, "Genus"]

norm_f.clr <- microbiome::transform(x = ps_its, transform = "clr") %>% 
    otu_table() %>% t() %>% data.frame()
colnames(norm_f.clr) <- data.frame( tax_table( object = ps_its ) )[, "Genus"]

norm_m.clr <- as.data.frame( compositions::clr(metabo+0.1))

# Correlations ------------------------------------------------------------

# 16S vs Metabolites
dat <- cbind(
    norm_b.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "Low fiber diet", ], 
    norm_m.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "Low fiber diet",])

corr_low <- cor(dat)
pval_low <- cor.mtest(mat = dat, method = "spearman")$p
corr_low <- corr_low[ colnames(norm_m.clr), colnames(norm_b.clr) ]
pval_low <- pval_low[ colnames(norm_m.clr), colnames(norm_b.clr)]

pval_tmp <- pval_low < 0.05
pval_tmp <- colSums(pval_tmp)
pval_tmp <- pval_tmp[pval_tmp > 0]

corr_low <- corr_low[ , names(pval_tmp) ]
pval_low <- pval_low[ , names(pval_tmp)]

dat <- cbind(
    norm_b.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "High fiber diet", ], 
    norm_m.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "High fiber diet",])

corr_high <- cor(dat)
pval_high <- cor.mtest(mat = dat, method = "spearman")$p
corr_high <- corr_high[ colnames(norm_m.clr), colnames(norm_b.clr) ]
pval_high <- pval_high[ colnames(norm_m.clr), colnames(norm_b.clr)]

pval_tmp <- pval_high < 0.05
pval_tmp <- colSums(pval_tmp)
pval_tmp <- pval_tmp[pval_tmp > 0]

corr_high <- corr_high[ , names(pval_tmp) ]
pval_high <- pval_high[ , names(pval_tmp)]

df_low <- corr_low %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, rho = value) %>% 
    dplyr::mutate(Group = "LFD")
tmp <- pval_low %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, pval = value) %>% 
    dplyr::mutate(Group = "LFD")
df_low$pvalue <- tmp$pval 

df_high <- corr_high %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, rho = value) %>% 
    dplyr::mutate(Group = "HFD")
tmp <- pval_high %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, pval = value) %>% 
    dplyr::mutate(Group = "HFD")
df_high$pvalue <- tmp$pval 

df_16S <- rbind(df_low, df_high)

# ITS vs Metabolites
dat <- cbind(norm_f.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "Low fiber diet",], 
             norm_m.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "Low fiber diet",])
corr_low <- cor(dat)
pval_low <- cor.mtest(mat = dat, method = "spearman")$p
corr_low <- corr_low[ colnames(norm_m.clr), colnames(norm_f.clr) ]
pval_low <- pval_low[ colnames(norm_m.clr), colnames(norm_f.clr)]

pval_tmp <- pval_low < 0.05
pval_tmp <- colSums(pval_tmp)
pval_tmp <- pval_tmp[pval_tmp > 0]

corr_low <- corr_low[ , names(pval_tmp) ]
pval_low <- pval_low[ , names(pval_tmp)]

dat <- cbind(norm_f.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "High fiber diet",], 
             norm_m.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "High fiber diet",])
corr_high <- cor(dat)
pval_high <- cor.mtest(mat = dat, method = "spearman")$p
corr_high <- corr_high[ colnames(norm_m.clr), colnames(norm_f.clr) ]
pval_high <- pval_high[ colnames(norm_m.clr), colnames(norm_f.clr)]

pval_tmp <- pval_high < 0.05
pval_tmp <- colSums(pval_tmp)
pval_tmp <- pval_tmp[pval_tmp > 0]

corr_high <- corr_high[ , names(pval_tmp) ]
pval_high <- pval_high[ , names(pval_tmp)]

df_low <- corr_low %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, rho = value) %>% 
    dplyr::mutate(Group = "LFD")
tmp <- pval_low %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, pval = value) %>% 
    dplyr::mutate(Group = "LFD")
df_low$pvalue <- tmp$pval 

df_high <- corr_high %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, rho = value) %>% 
    dplyr::mutate(Group = "HFD")
tmp <- pval_high %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Metabolite = Var1, Genus = Var2, pval = value) %>% 
    dplyr::mutate(Group = "HFD")
df_high$pvalue <- tmp$pval 

df_ITS <- rbind(df_low, df_high)

# Correlations 16S vs ITS -------------------------------------------------

dat <- cbind(norm_b.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "Low fiber diet",], 
             norm_f.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "Low fiber diet",])
corr_low <- cor(dat)
pval_low <- cor.mtest(mat = dat, method = "spearman")$p
corr_low <- corr_low[ colnames(norm_f.clr), colnames(norm_b.clr) ]
pval_low <- pval_low[ colnames(norm_f.clr), colnames(norm_b.clr)]

pval_tmp <- pval_low < 0.05
pval_tmp <- colSums(pval_tmp)
pval_tmp <- pval_tmp[pval_tmp > 0]

corr_low <- corr_low[ , names(pval_tmp) ]
pval_low <- pval_low[ , names(pval_tmp)]

dat <- cbind(norm_b.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "High fiber diet",], 
             norm_f.clr[cova$Status=="control" & cova$Fiber.Content.Diet == "High fiber diet",])
corr_high <- cor(dat)
pval_high <- cor.mtest(mat = dat, method = "spearman")$p
corr_high <- corr_high[ colnames(norm_f.clr), colnames(norm_b.clr) ]
pval_high <- pval_high[ colnames(norm_f.clr), colnames(norm_b.clr)]

pval_tmp <- pval_high < 0.05
pval_tmp <- colSums(pval_tmp)
pval_tmp <- pval_tmp[pval_tmp > 0]

corr_high <- corr_high[ , names(pval_tmp) ]
pval_high <- pval_high[ , names(pval_tmp)]

df_low <- corr_low %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Fungi = Var1, Genus = Var2, rho = value) %>% 
    dplyr::mutate(Group = "LFD")
tmp <- pval_low %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Fungi = Var1, Genus = Var2, pval = value) %>% 
    dplyr::mutate(Group = "LFD")
df_low$pvalue <- tmp$pval 

df_high <- corr_high %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Fungi = Var1, Genus = Var2, rho = value) %>% 
    dplyr::mutate(Group = "HFD")
tmp <- pval_high %>% 
    reshape2::melt(data = .) %>% 
    dplyr::select(Fungi = Var1, Genus = Var2, pval = value) %>% 
    dplyr::mutate(Group = "HFD")
df_high$pvalue <- tmp$pval 

df_16S_ITS <- rbind(df_low, df_high)

# Plots -------------------------------------------------------------------

p_corr_bact <- df_16S %>% 
    dplyr::filter( abs(rho) > 0.3 ) %>%
    dplyr::filter( pvalue < 0.01 ) %>% 
    arrange(Metabolite) %>% 
    dplyr::mutate(Metabolite = gsub(pattern = "U6\\.", replacement = "U6", x = Metabolite)) %>% 
    dplyr::mutate(Metabolite = gsub(pattern = "\\.", replacement = " ", x = Metabolite)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "\\.", replacement = " ", x = Genus)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus)) %>% 
    ggplot(data = ., aes(x = Genus, 
                         y =factor(Metabolite, levels = rev(levels(factor(Metabolite)))), 
                         fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, 
                         limit = c(-1,1), space = "Lab",
                         name="Spearman's rho") +
    labs(title = "Correlation Matrix 16S", 
         x = "", y = "") +
    facet_grid(facets = ~ Group,  
               scales = "free_x", space = "free_x") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=8, face = "italic"),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_corr_bact

p_corr_fung <- df_ITS %>% 
    dplyr::filter( abs(rho) > 0.3 ) %>%
    dplyr::filter( pvalue < 0.01 ) %>% 
    arrange(Metabolite) %>% 
    dplyr::mutate(Metabolite = gsub(pattern = "U6\\.", replacement = "U6", x = Metabolite)) %>% 
    dplyr::mutate(Metabolite = gsub(pattern = "\\.", replacement = " ", x = Metabolite)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "\\.", replacement = " ", x = Genus)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus)) %>% 
    ggplot(data = ., aes(x = Genus, 
                         y =factor(Metabolite, levels = rev(levels(factor(Metabolite)))), 
                         fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, 
                         limit = c(-1,1), space = "Lab",
                         name="Spearman's rho") +
    labs(title = "Correlation Matrix ITS", 
         x = "", y = "") +
    facet_grid(facets = ~ Group,  
               scales = "free_x", space = "free_x") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=8, face = "italic"),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          # axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_corr_fung

p_corr_bact_fung <- df_16S_ITS %>% 
    dplyr::filter( abs(rho) > 0.3 ) %>%
    dplyr::filter( pvalue < 0.01 ) %>% 
    arrange(Fungi) %>% 
    dplyr::mutate(Fungi = gsub(pattern = "_", replacement = " ", x = Fungi)) %>% 
    dplyr::mutate(Fungi = gsub(pattern = "\\.", replacement = " ", x = Fungi)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "\\.", replacement = " ", x = Genus)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus)) %>% 
    ggplot(data = ., aes(x = Genus, 
                         y =factor(Fungi, levels = rev(levels(factor(Fungi)))), 
                         fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, 
                         limit = c(-1,1), space = "Lab",
                         name="Spearman's rho") +
    labs(title = "Correlation Matrix 16S vs ITS", 
         x = "", y = "") +
    facet_grid(facets = ~ Group,  
               scales = "free_x", space = "free_x") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=8, face = "italic"),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=8, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_corr_bact_fung + ggtitle('Correlations 16S vs ITS ( p < 0.01)')

layout <- "
AAACCCCC
BBXCCCCC
"
p1 <- p_corr_bact + theme(legend.position = "none") + labs(title = "Bacteriome vs Metabolites", x = "", y = "") +
    p_corr_fung + theme(legend.position = "none") + labs(title = "Mycobiome vs Metabolites", x = "", y = "") +
    p_corr_bact_fung + theme(legend.position = "right") + labs(title = "Bacteriome vs Mycobiome", x = "", y = "") +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Correlations healthy cohort (p < 0.01)') 
p1 + plot_annotation(title = '') 
ggsave(filename = "plots_final/Suppl_Fig_6.pdf", width = 27, height = 9, 
       plot = p1 + plot_annotation(title = '') )
