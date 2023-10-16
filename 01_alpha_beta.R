# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)

library(ggpubr)
library(patchwork)

# Colors ------------------------------------------------------------------

# Colors for taxonomic abundance plots
taxNumberPlot <- 19
cols <- c(pals::tableau20(taxNumberPlot+1), "grey85")

colv <- c("steelblue", "orange")

# Get processed data ------------------------------------------------------

ps_16s <- readRDS(file = "processed_ps_16s.RDS")
ps_its <- readRDS(file = "processed_ps_its.RDS")
metabo <- readRDS(file = "processed_metabo.RDS")
cova   <- readRDS(file = "processed_cova.RDS")

suppl_table <- cova %>% 
    dplyr::select(Sample, Status, Gender, Age, 
                  Diet_Fiber = Fiber.Content.Diet, BMI = BMI..Kg.m) %>% 
    dplyr::mutate(Status = case_when(
        Status == 'control' ~ 'Healthy',
        Status == 'disease' ~ 'T2D',
    )) %>% 
    dplyr::mutate(Diet_Fiber = case_when(
        Diet_Fiber == 'High fiber diet' ~ 'High',
        Diet_Fiber == 'Low fiber diet' ~ 'Low',
    ))

WriteXLS::WriteXLS(x = suppl_table, ExcelFileName = "Suppl_Table_cova.xlsx", 
                   SheetNames = 'Cohort', 
                   row.names = F, AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)

ps_met <-phyloseq(otu_table(object = metabo %>% 
                                column_to_rownames("Sample"), 
                            taxa_are_rows = F), 
                  sample_data(ps_16s) )

sample_data(ps_met)$Shannon_Metabolome <- estimate_richness(physeq = ps_met, measures = "Shannon")$Shannon

sample_data(ps_16s)$Disease <- factor( gsub(pattern = "Healthy", replacement = "Control", x = sample_data(ps_16s)$Disease) )
sample_data(ps_its)$Disease <- factor( gsub(pattern = "Healthy", replacement = "Control", x = sample_data(ps_its)$Disease) )
sample_data(ps_met)$Disease <- factor( gsub(pattern = "Healthy", replacement = "Control", x = sample_data(ps_met)$Disease) )

sample_data(ps_16s)$Var <- paste(sample_data(ps_16s)$Disease, sample_data(ps_16s)$Diet)
sample_data(ps_its)$Var <- paste(sample_data(ps_its)$Disease, sample_data(ps_its)$Diet)
sample_data(ps_met)$Var <- paste(sample_data(ps_met)$Disease, sample_data(ps_met)$Diet)

# Age ---------------------------------------------------------------------

table(cova$Status, cova$Fiber.Content.Diet)
table(cova$Gender, cova$Fiber.Content.Diet)

p_age_1 <- ps_16s %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., x = "Disease", y = "Age", 
                     trim = TRUE, width = 0.25,
                     palette = c("grey35", "grey65"), fill = "Disease", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    ggpubr::stat_compare_means() +
    ylim(0,100) +
    xlab("") + ylab("Age") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_age_1

p_age_2 <- ps_16s %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., x = "Diet", y = "Age", 
                     trim = TRUE, width = 0.25,
                     palette = colv, fill = "Diet", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    facet_wrap(Disease~., scales = "free_x") +
    ggpubr::stat_compare_means() +
    ylim(0,100) +
    xlab("") + ylab("Age") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_age_2
ggsave(filename = "plots/age.pdf", width = 4, height = 6)


p_bmi_1 <- ps_16s %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., x = "Disease", y = "BMI", 
                     trim = TRUE, width = 0.25,
                     palette = colv, fill = "Disease", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    ggpubr::stat_compare_means() +
    ylim(0,75) +
    xlab("") + ylab("BMI") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_bmi_1

p_bmi_2 <- ps_16s %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., x = "Diet", y = "BMI", 
                     trim = TRUE, width = 0.25,
                     palette = colv, fill = "Diet", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    facet_wrap(Disease~., scales = "free_x") +
    ggpubr::stat_compare_means() +
    ylim(0,75) +
    xlab("") + ylab("BMI") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_bmi_2

ps_16s %>% 
    sample_data() %>% 
    data.frame() %>% 
    group_by(Disease) %>% 
    summarise(mean_BMI = mean(BMI), sd_BMI = sd(BMI)) 
# Top 3 taxa --------------------------------------------------------------

ps_16s %>%
    microbiome::aggregate_taxa(x = ., level = "Genus") %>% 
    otu_table() %>% t() %>%
    colSums() %>% sort(decreasing = T) %>% head( n = 5)

ps_its %>%
    microbiome::aggregate_taxa(x = ., level = "Genus") %>% 
    otu_table() %>% t() %>%
    colSums() %>% sort(decreasing = T) %>% head( n = 5)

ps_met %>%
    otu_table() %>% 
    colSums() %>% sort(decreasing = T) %>% head( n = 5)

# 16S Taxa Plots ----------------------------------------------------------

p16s_phylum <- ps_16s %>%
    tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>%                                         # Melt to long format
    #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
    arrange(Phylum) %>%
    dplyr::select( Phylum, Diet, Disease, Abundance) %>%
    dplyr::group_by(Phylum, Diet, Disease) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Phylum)) %>%
    unique() %>%
    dplyr::mutate(Phylum = gsub(pattern = "p__", replacement = "", x = Phylum)) %>%
    ggplot(., aes(x = Diet, y = mean, fill = Phylum)) +
    facet_wrap(Disease~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance \n") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
p16s_phylum

# Genera
means_genus <- ps_16s %>%
    tax_glom(taxrank = "Genus") %>% # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>% # Melt to long format
    dplyr::arrange(Genus) %>%
    dplyr::select( Genus, Diet, Disease, Abundance) %>%
    dplyr::group_by(Genus, Diet, Disease) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Genus)) %>%
    unique() %>%
    dplyr::mutate(Genus = gsub(pattern = "g__", replacement = "", x = Genus))
# filter(mean >= 0.02) # Filter out low abundance taxa

# summarise lowly abundant genera
means_genus$Genus[means_genus$mean <= 0.02] <- "Other"
means_genus <- means_genus %>%
    group_by(Genus, Diet, Disease) %>%
    summarise(mean = sum(mean), groups = c(Genus)) %>%
    unique()

means_genus$mean[means_genus$Disease=="DM II" & means_genus$Diet == "LFD" & means_genus$Genus == "Other"] <- 
    means_genus$mean[means_genus$Disease=="DM II" & means_genus$Diet == "LFD" & means_genus$Genus == "Other"] + 
    ( 1 - sum(means_genus$mean[means_genus$Disease=="DM II" & means_genus$Diet == "LFD"]) )
means_genus$mean[means_genus$Disease=="DM II" & means_genus$Diet == "HFD" & means_genus$Genus == "Other"] <- 
    means_genus$mean[means_genus$Disease=="DM II" & means_genus$Diet == "HFD" & means_genus$Genus == "Other"] + 
    ( 1 - sum(means_genus$mean[means_genus$Disease=="DM II" & means_genus$Diet == "HFD"]) )

means_genus$mean[means_genus$Disease=="Healthy" & means_genus$Diet == "LFD" & means_genus$Genus == "Other"] <- 
    means_genus$mean[means_genus$Disease=="Healthy" & means_genus$Diet == "LFD" & means_genus$Genus == "Other"] + 
    ( 1 - sum(means_genus$mean[means_genus$Disease=="Healthy" & means_genus$Diet == "LFD"]) )
means_genus$mean[means_genus$Disease=="Healthy" & means_genus$Diet == "HFD" & means_genus$Genus == "Other"] <- 
    means_genus$mean[means_genus$Disease=="Healthy" & means_genus$Diet == "HFD" & means_genus$Genus == "Other"] + 
    ( 1 - sum(means_genus$mean[means_genus$Disease=="Healthy" & means_genus$Diet == "HFD"]) )

p16s_genus <- means_genus %>% 
    dplyr::mutate(Genus = gsub(pattern = "f__", replacement = "", x = Genus)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "_unclassified", replacement = " uncl.", x = Genus)) %>% 
    ggplot(data = ., aes(x = Diet, y = mean, fill = Genus)) +
    facet_wrap(Disease~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance (Genus abundance > 2%)") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Genus", ncol = 1))
p16s_genus

# 16S Alpha ---------------------------------------------------------------

tmp <- sample_data(ps_16s) %>% data.frame() 
m1 <- glm(Shannon_16s ~ Age, data = tmp)
tmp$Residuals <- residuals(m1)
fitted(m1)
residuals(m1)

p16s_alpha <- ps_16s %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., 
                     x = "Diet", y = "Shannon_16s", 
                     trim = TRUE, width = 0.5,
                     palette = colv, fill = "Diet", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    facet_wrap(Disease~., scales = "free_x") +
    ggpubr::stat_compare_means() +
    ylim(0,6) +
    xlab("") + ylab("Shannon index H") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p16s_alpha

# 16S Beta diversity ------------------------------------------------------

ps_genus <- tax_glom(physeq = ps_16s, taxrank = "Genus")
ps_ra <- transform_sample_counts(ps_genus, function(x){x / sum(x)})
ps_genus_clr <- microbiome::transform(ps_genus, "clr")

sample_data(ps_genus_clr)$Prevotella <- ps_ra %>% 
    microbiome::aggregate_taxa(x = ., level = "Genus") %>%
    subset_taxa(physeq = ., Genus == "g__Prevotella") %>% 
    otu_table() %>% 
    data.frame() %>% t() 

ps_clr <- microbiome::transform(ps_16s, "clr")
sample_data(ps_clr)$Prevotella <- sample_data(ps_genus_clr)$Prevotella
otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")
beta_16S <- ps_clr_dist
ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
p16s_beta <- plot_ordination( physeq = ps_clr, ordination = ps_clr_ord, 
                              color='Prevotella', shape = "Var")  +
    # scale_color_manual(values=colv) +
    scale_color_continuous(type = "viridis") +
    geom_point(size=3) +
    ggtitle( paste0("RDA (Aitchison distance; genus level)" ) ) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p16s_beta

beta_df <- vegan::adonis2(formula = ps_clr_dist ~ Diet * Disease + Prevotella + Age + Gender + BMI,
                         data = data.frame(ps_clr %>% sample_data),
                         permutations = 99999) %>% data.frame() %>%
    rownames_to_column("Variable")
beta_df

# Otu00001 -> Prevotella (tax_table(ps_genus))
# Null model
corncob_null <- corncob::bbdml(formula = Otu00001 ~ Shannon_16s, 
                             phi.formula = ~ Shannon_16s,
                             data = ps_genus %>% subset_samples(physeq = ., Disease == "Control"))
# DA Model
corncob_da <- corncob::bbdml(formula = Otu00001 ~ Diet + Shannon_16s, 
                          phi.formula = ~ Diet + Shannon_16s,
                          data = ps_genus %>% subset_samples(physeq = ., Disease == "Control"))
plot(corncob_da, color = "Var", total = F, B = 50)

corncob::lrtest(mod_null = corncob_null, mod = corncob_da)
corncob::waldchisq(mod = corncob_da, mod_null = corncob_null) 
corncob::waldt(mod = corncob_da) %>% data.frame()

# Null model
corncob_null <- corncob::bbdml(formula = Otu00001 ~ Shannon_16s, 
                               phi.formula = ~ Shannon_16s,
                               data = ps_genus %>% subset_samples(physeq = ., Disease == "DM II"))
# DA Model
corncob_da <- corncob::bbdml(formula = Otu00001 ~ Diet + Shannon_16s, 
                             phi.formula = ~ Diet + Shannon_16s,
                             data = ps_genus %>% subset_samples(physeq = ., Disease == "DM II"))
plot(corncob_da, color = "Diet", total = F, B = 50)

corncob::lrtest(mod_null = corncob_null, mod = corncob_da)
corncob::waldchisq(mod = corncob_da, mod_null = corncob_null) 

df_waldt <- corncob::waldt(mod = corncob_da) %>% data.frame()
df_waldt

p16s_prevot <- p16s_beta$data %>% 
    data.frame() %>% 
    ggboxplot(data = ., x = "Diet", y = "Prevotella", 
              outlier.shape = NA, alpha = 0.75,
              add = "none", width = 0.3, fill = "Diet") + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    scale_fill_manual(values=colv) +
    facet_wrap(Disease~., scales = "free_x") +
    # scale_y_continuous(position = "right", sec.axis = sec_axis(~., labels = NULL)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    labs(x = "", y = "Rel. abundances Prevotella") + ylim(0,1) 
p16s_prevot

cor.test(p16s_beta$data$PC1, p16s_beta$data$Shannon_16s, method = 'spearman')
cor.test(p16s_beta$data$PC2, p16s_beta$data$Shannon_16s, method = 'spearman')
cor.test(p16s_beta$data$PC1, p16s_beta$data$Prevotella, method = 'spearman')
cor.test(p16s_beta$data$PC2, p16s_beta$data$Prevotella, method = 'spearman')

m1 <- lm(PC1 ~ Diet * Disease + Prevotella + Shannon_16s + Shannon_ITS + Age + Gender + BMI, data = p16s_beta$data)
m2 <- lm(PC2 ~ Diet * Disease + Prevotella + Shannon_16s + Shannon_ITS + Age + Gender + BMI, data = p16s_beta$data)
summary(m1)
summary(m2)

wilcox.test(Prevotella ~ Disease, data = p16s_beta$data)

layout <- "
ABCCDD
ABEEEE
"

p16s_phylum + p16s_genus + 
    p16s_alpha + p16s_prevot +
    p16s_beta +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = '16S rRNA gene data') 
ggsave(filename = "plots/plots_16s.pdf", height = 9, width = 18)

# ITS Taxa plots ----------------------------------------------------------

pits_phylum <- ps_its %>%
    tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>%                                         # Melt to long format
    #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
    arrange(Phylum) %>%
    dplyr::select( Phylum, Diet, Disease, Abundance) %>%
    dplyr::group_by(Phylum, Diet, Disease) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Phylum)) %>%
    unique() %>%
    dplyr::mutate(Phylum = gsub(pattern = "p__", replacement = "", x = Phylum)) %>%
    ggplot(., aes(x = Diet, y = mean, fill = Phylum)) +
    facet_wrap(Disease~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance \n") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
pits_phylum

# Genera
means_genus <- ps_its %>%
    tax_glom(taxrank = "Genus") %>% # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>% # Melt to long format
    dplyr::arrange(Genus) %>%
    dplyr::select( Genus, Diet, Disease, Abundance) %>%
    dplyr::group_by(Genus, Diet, Disease) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Genus)) %>%
    unique() %>%
    dplyr::mutate(Genus = gsub(pattern = "g__", replacement = "", x = Genus))
# filter(mean >= 0.02) # Filter out low abundance taxa

# summarise lowly abundant genera
means_genus$Genus[means_genus$mean <= 0.01] <- "Other"
means_genus <- means_genus %>%
    group_by(Genus, Diet, Disease) %>%
    summarise(mean = sum(mean), groups = c(Genus)) %>%
    unique()

# means_genus$mean[means_genus$Disease == "DM II" & means_genus$Diet == "HFD"] %>% sum()

pits_genus <- means_genus %>% 
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus)) %>% 
    dplyr::mutate(Genus = gsub(pattern = "unclassified", replacement = "uncl.", x = Genus)) %>% 
    ggplot(data = ., aes(x = Diet, y = mean, fill = Genus)) +
    facet_wrap(Disease~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance (Genus abundance > 1%)") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Genus", ncol = 1))
pits_genus

# ITS Alpha ---------------------------------------------------------------

tmp <- sample_data(ps_its) %>% data.frame() 
m1 <- glm(Shannon_ITS ~ Age, data = tmp)
tmp$Residuals <- residuals(m1)
fitted(m1)
residuals(m1)

pits_alpha <- ps_its %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., 
                     x = "Diet", y = "Shannon_ITS", 
                     trim = TRUE, width = 0.5,
                     palette = colv, fill = "Diet", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    facet_wrap(Disease~., scales = "free_x") +
    ggpubr::stat_compare_means() +
    ylim(0,6) +
    xlab("") + ylab("Shannon index H") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
pits_alpha

# ITS Beta diversity ------------------------------------------------------

# ps_genus <- tax_glom(physeq = ps_its, taxrank = "Genus")
ps_clr <- microbiome::transform(ps_its, "clr")

sample_data(ps_clr)$Prevotella <- ps_ra %>% 
    microbiome::aggregate_taxa(x = ., level = "Genus") %>%
    subset_taxa(physeq = ., Genus == "g__Prevotella") %>% 
    otu_table() %>% 
    data.frame() %>% t() 

otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")
beta_ITS <- ps_clr_dist
ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
pits_beta <- plot_ordination( physeq = ps_clr, ordination = ps_clr_ord, 
                              color='Diet', shape = "Disease")  +
    scale_color_manual(values=colv) +
    # scale_color_continuous(type = "viridis") +
    geom_point(size=3) +
    ggtitle( paste0("RDA of Aitchison distance (Genus level)" ) ) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
pits_beta

beta_df_its <- vegan::adonis2(formula = ps_clr_dist ~ Diet * Disease + Prevotella + Age +Gender + BMI,
                             data = data.frame(ps_clr %>% sample_data),
                             permutations = 99999) %>% data.frame() %>%
    rownames_to_column("Variable")
beta_df_its

m1 <- lm(PC1 ~ Diet * Disease + Prevotella + Shannon_ITS + Shannon_16s + Age + Gender + BMI, data = pits_beta$data)
m2 <- lm(PC2 ~ Diet * Disease + Prevotella + Shannon_ITS + Shannon_16s + Age + Gender + BMI, data = pits_beta$data)
summary(m1)
summary(m2)

layout <- "
ABCCCC
ABDDDD
"
pits_phylum + pits_genus + 
    pits_alpha + pits_beta +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'ITS data') 
ggsave(filename = "plots/plots_its.pdf", height = 9, width = 18)

# Alpha 16S vs Alpha ITS --------------------------------------------------

alpha <- sample_data(ps_16s) %>% data.frame() 

cor.test(alpha$Shannon_16s, alpha$Shannon_ITS)

summary( lm(Shannon_16s~Shannon_ITS, data = alpha))

alpha %>% 
    ggplot(.,aes(x = Shannon_16s, y = Shannon_ITS)) +
    geom_point() +
    geom_smooth(method='lm') +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    xlim(2,5) + ylim(0,3)

# Metabolome abundance plots --------------------------------------------------

# Abundances
means_metab <- ps_met %>%
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>% # Melt to long format
    dplyr::rename(Metabolite = OTU) %>% 
    dplyr::arrange(Metabolite) %>%
    dplyr::select( Metabolite, Diet, Disease, Abundance) %>%
    dplyr::group_by(Metabolite, Diet, Disease) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Metabolite)) %>%
    unique() %>%
    dplyr::mutate(Metabolite = gsub(pattern = "g__", replacement = "", x = Metabolite))
# filter(mean >= 0.02) # Filter out low abundance taxa

# summarise lowly abundant genera
means_metab$Metabolite[means_metab$mean < 0.02] <- "Other"
means_metab <- means_metab %>%
    group_by(Metabolite, Diet, Disease) %>%
    summarise(mean = sum(mean), groups = c(Metabolite)) %>%
    unique()

pmet_abund <- means_metab %>% 
    dplyr::mutate(Metabolite = gsub(pattern = "\\.", replacement = "", x = Metabolite)) %>%
    dplyr::mutate(Metabolite = gsub(pattern = "X2", replacement = "2-", x = Metabolite)) %>%
    ggplot(data = ., aes(x = Diet, y = mean, fill = Metabolite)) +
    facet_wrap(Disease~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance (Metabolite abundance > 2%)") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Metabolite", ncol = 1))
pmet_abund

# Metabolome Alpha ------------------------------------------------------------

pmet_alpha <- ps_met %>% 
    sample_data() %>% 
    data.frame() %>% 
    ggpubr::ggviolin(data = ., 
                     x = "Diet", y = "Shannon_Metabolome", 
                     trim = TRUE, width = 0.5,
                     palette = colv, fill = "Diet", 
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(color="grey25")) +
    facet_wrap(Disease~., scales = "free_x") +
    ggpubr::stat_compare_means() +
    ylim(0,6) +
    xlab("") + ylab("Shannon index H") + 
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
pmet_alpha

# Metabolome Beta diversity -------------------------------------------------

ps_clr <- microbiome::transform(ps_met, "clr")

sample_data(ps_clr)$Prevotella <- ps_ra %>% 
    microbiome::aggregate_taxa(x = ., level = "Genus") %>%
    subset_taxa(physeq = ., Genus == "g__Prevotella") %>% 
    otu_table() %>% 
    data.frame() %>% t() 

otu.table_clr <- otu_table(ps_clr) 
ps_clr_dist <- dist(otu.table_clr, method="euclidean")
beta_Met <- ps_clr_dist
ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
pmet_beta <- plot_ordination( physeq = ps_clr, ordination = ps_clr_ord, 
                              color='Diet', shape = "Disease")  +
    scale_color_manual(values=colv) +
    # scale_color_continuous(type = "viridis") +
    geom_point(size=3) +
    ggtitle( paste0("RDA of Aitchison distance" ) ) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(),
          legend.key = element_rect(fill = "transparent"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
pmet_beta

beta_df_met <- vegan::adonis2(formula = ps_clr_dist ~ Diet * Disease + Prevotella + Age + Gender + BMI,
                              data = data.frame(ps_clr %>% sample_data),
                              permutations = 99999) %>% data.frame() %>%
    rownames_to_column("Variable")
beta_df_met

m1 <- lm(PC1 ~ Diet * Disease + Prevotella + Shannon_Metabolome + Age + Gender + BMI, data = pmet_beta$data)
m2 <- lm(PC2 ~ Diet * Disease + Prevotella + Shannon_Metabolome + Age + Gender + BMI, data = pmet_beta$data)
summary(m1)
summary(m2)

layout <- "
ABB
ACC
"
pmet_abund + 
    pmet_alpha + pmet_beta +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Metabolome data') 
ggsave(filename = "plots/plots_metabolome.pdf", height = 9, width = 9)

# Compare diversity -------------------------------------------------------

vegan::mantel(xdis = beta_16S, ydis = beta_ITS, method = "spearman", 
              permutations = 10000, parallel = 4)
# Mantel statistic r: -0.1203 
# Significance: 0.89241 
vegan::mantel(xdis = beta_16S, ydis = beta_Met, method = "spearman", 
              permutations = 10000, parallel = 4)
# Mantel statistic r: 0.2095 
# Significance: 0.012599 
vegan::mantel(xdis = beta_ITS, ydis = beta_Met, method = "spearman", 
              permutations = 10000, parallel = 4)
# Mantel statistic r: -0.03797 
# Significance: 0.62504 

cor.test(x = sample_data(ps_16s)$Shannon_16s, 
         y = sample_data(ps_its)$Shannon_ITS, 
         method = "spearman")

cor.test(x = sample_data(ps_16s)$Shannon_16s, 
         y = sample_data(ps_met)$Shannon_Met, 
         method = "spearman")

cor.test(x = sample_data(ps_met)$Shannon_Met, 
         y = sample_data(ps_its)$Shannon_ITS, 
         method = "spearman")

# Pretty plotting ---------------------------------------------------------

# Figure 1
layout <- "
ABC
DEF
"
p16s_alpha + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    pits_alpha + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    pmet_alpha + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    p16s_beta + ggtitle("") + theme(legend.justification = "top") + 
        guides(
        colour = guide_legend(title = "Prevotella", title.theme = element_text(face = "italic")), 
        shape = guide_legend(title = "")) +
    pits_beta + ggtitle("") + theme(legend.position = "none") +
    pmet_beta + ggtitle("") + theme(legend.position = "none") +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Figure_1.pdf", height = 9, width = 18)

# Figure 2
layout <- "
ABC
"
p16s_genus + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    pits_genus + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    pmet_abund + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Figure_2.pdf", height = 7.5, width = 15)

# Suppl Fig 1
layout <- "
ABB
"
p_age_1 + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    p_age_2 + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Suppl_Fig_1.pdf", height = 4.5, width = 9)

# Suppl Fig 2
p16s_prevot + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) 
ggsave(filename = "plots_final/Suppl_Fig_2.pdf", height = 4.5, width = 4.5)

# Figure 3
layout <- "
AB
"
p16s_phylum + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    pits_phylum + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Suppl_Fig_3.pdf", height = 7.5, width = 10)

layout <- "
ABB
"
p_bmi_1 + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    p_bmi_2 + theme( axis.text.x = element_text(angle = 0, hjust = 0.5, size=12) ) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Suppl_Fig_BMI.pdf", height = 4.5, width = 9)
