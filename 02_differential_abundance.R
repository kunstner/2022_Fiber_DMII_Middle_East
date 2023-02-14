# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)

library(ggpubr)
library(patchwork)

library(ANCOMBC)
library(ALDEx2)
library(Maaslin2)

# devtools::install_github("barakbri/dacomp", force = T, build_vignettes = TRUE)
library(dacomp)
# vignette('dacomp_main_vignette')

# Colors ------------------------------------------------------------------

# Colors for taxonomic abundance plots
taxNumberPlot <- 19
cols <- c(pals::tableau20(taxNumberPlot+1), "grey85")

colv <- c("steelblue", "orange")

qcutoff <- 0.1

# Get processed data ------------------------------------------------------

ps_16s <- readRDS(file = "data/processed_ps_16s.RDS")
ps_its <- readRDS(file = "data/processed_ps_its.RDS")
metabo <- readRDS(file = "data/processed_metabo.RDS")
cova   <- readRDS(file = "data/processed_cova.RDS")

ps_met <-phyloseq(otu_table(object = metabo %>% 
                                column_to_rownames("Sample"), 
                            taxa_are_rows = F), 
                  sample_data(ps_16s) )

sample_data(ps_met)$Shannon_Metabolome <- estimate_richness(physeq = ps_met, measures = "Shannon")$Shannon

sample_data(ps_16s)$Var <- paste(sample_data(ps_16s)$Disease, sample_data(ps_16s)$Diet)
sample_data(ps_its)$Var <- paste(sample_data(ps_its)$Disease, sample_data(ps_its)$Diet)
sample_data(ps_met)$Var <- paste(sample_data(ps_met)$Disease, sample_data(ps_met)$Diet)

# Process data on Genus level ---------------------------------------------

# ps_16s_g <- tax_glom(physeq = ps_16s, taxrank = "Genus")
# ps_its_g <- tax_glom(physeq = ps_its, taxrank = "Genus")

# preserve genus names
ps_16s_g <- microbiome::aggregate_taxa(x = ps_16s, level = "Genus")
ps_its_g <- microbiome::aggregate_taxa(x = ps_its, level = "Genus")

# Prevalence filtering (20%)
ps_16s_g <- metagMisc::phyloseq_filter_prevalence(physeq = ps_16s_g, prev.trh = 0.2)
ps_its_g <- metagMisc::phyloseq_filter_prevalence(physeq = ps_its_g, prev.trh = 0.2)
ps_met_g <- metagMisc::phyloseq_filter_prevalence(physeq = ps_met, prev.trh = 0.2)

# 16S: DA Control Group --------------------------------------------------------

ps_16s_data <- ps_16s_g %>% subset_samples(physeq = ., Disease == "Healthy")

#
# ALDEx2
#

da_aldex <- aldex.clr(
    reads = ps_16s_data %>% otu_table(object = ., taxa_are_rows = T) %>% data.frame(),
    conds = ps_16s_data %>% sample_data(object = .) %>% .$Diet, 
    mc.samples = 512, 
    denom = "all",
    verbose = TRUE
)

# calculate expected values of the Welch's t-test and Wilcoxon rank test 
da_aldex_tt <- aldex.ttest(
    clr = da_aldex, 
    paired.test = FALSE, 
    verbose = TRUE)
# effect sizes
da_aldex_effect <- aldex.effect(da_aldex, CI = TRUE, verbose = FALSE)
# combine outputs 
aldex_16s_healthy <- data.frame(da_aldex_tt, da_aldex_effect)

aldex_16s_healthy %>% 
    rownames_to_column("genus") %>%
    dplyr::filter(we.eBH <= qcutoff)  %>% # t-test output
    dplyr::select(genus, we.eBH, we.eBH, effect, overlap) %>% 
    knitr::kable()

#
# ANCOMBC
#

# perform the analysis 
da_ancombc <- ancombc(
    phyloseq = ps_16s_data, 
    formula = "Diet", group = "Diet", tax_level = "Genus",
    p_adj_method = "BH", alpha = qcutoff, 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = TRUE, # recommended if the sample size is small 
    global = FALSE
)
# store the results in res 
da_ancombc_healthy <- da_ancombc$res

da_ancombc_healthy$diff_abn %>% 
    dplyr::filter(DietLFD == TRUE) %>% 
    knitr::kable()

#
# MaAsLin2
#

# similar settings as in Nearing et al. (2021)
da_maaslin <- Maaslin2( 
    input_data = ps_16s_data %>% otu_table(object = ., taxa_are_rows = T) %>% 
        data.frame() %>% t(), 
    input_metadata = ps_16s_data %>% sample_data() %>% data.frame(),
    output = "MaAsLin2_16s_healthy",
    transform = "AST", # AST: arcsine square root
    fixed_effects = "Diet",
    reference = "Diet,HFD",  
    normalization = "TSS",
    standardize = FALSE, cores = 4, 
    plot_heatmap = T, heatmap_first_n = 10,
    min_prevalence = 0, # prev filterin already done
    correction = "BH"
)

da_maaslin$results %>% 
    dplyr::filter(qval < qcutoff)

#
# dacomp
#

df_dacomp <- ps_16s_data %>% 
    otu_table(object = ., taxa_are_rows = T) %>% 
    data.frame() %>% t()
result.selected.references <- dacomp.select_references(
    X = df_dacomp,
    minimal_TA = 100, 
    verbose = T)
length(result.selected.references$selected_references)
dacomp.plot_reference_scores(result.selected.references)
colnames(df_dacomp)[result.selected.references$selected_references]

#multiplicity correction levels for the BH and DS-FDR methods
q_BH <- q_DSFDR <- qcutoff

#Perform testing:
result.test <- dacomp.test(
    X = df_dacomp,
    y = ps_16s_data %>% 
        sample_data(object = .) %>% .$Diet,
    ind_reference_taxa = result.selected.references,
    test = DACOMP.TEST.NAME.WILCOXON,
    # nr_perm = 1000, # use default values here
    verbose = T, 
    q = q_DSFDR,
    return_rarefied_values = TRUE)

rejected_BH <- which(p.adjust(result.test$p.values.test,method = 'BH') <= q_BH)
rejected_DSFDR <- result.test$dsfdr_rejected

da_dacomp <- data.frame(
    feature = df_dacomp %>% colnames(),
    pval = result.test$p.values.test,
    qval_DSFDR = result.test$p.values.test.adjusted,
    effect = result.test$effect_size_estimates
)

#
# Summarise methods
#

da_16s_healthy <- dplyr::full_join(
    x = rownames_to_column(aldex_16s_healthy, "genus") %>%
        dplyr::select(genus, aldex2 = we.eBH),
    y = rownames_to_column(da_ancombc$res$diff_abn, "genus") %>%
        dplyr::select(genus, ancombc = DietLFD),
    by = "genus") %>%
    dplyr::full_join( x = .,
        y = dplyr::select(da_maaslin$results, genus = feature, maaslin2 = qval), 
        by = "genus") %>%
    dplyr::full_join( x = .,
        dplyr::select(da_dacomp, genus = feature, dacomp = qval_DSFDR), 
        by = "genus") %>%
    dplyr::mutate(
        across(c(aldex2, maaslin2, dacomp), ~ .x <= qcutoff),
        score = rowSums(across(c(aldex2, ancombc, maaslin2, dacomp)), na.rm = T)
    ) %>% 
    dplyr::filter( !is.na(score) & score > 0)
da_16s_healthy$aldex2[ is.na(da_16s_healthy$aldex2) ] <- FALSE
da_16s_healthy$ancombc[ is.na(da_16s_healthy$ancombc) ] <- FALSE
da_16s_healthy$maaslin2[ is.na(da_16s_healthy$maaslin2) ] <- FALSE
da_16s_healthy$dacomp[ is.na(da_16s_healthy$dacomp) ] <- FALSE

# how many genera were identified by each method?
summarise(da_16s_healthy, across(where(is.logical), sum)) %>%
    knitr::kable()

da_16s_healthy_2 <- da_16s_healthy %>% 
    dplyr::filter(score >= 2) %>% 
    dplyr::mutate(genus = gsub(pattern = "\\[|\\]", replacement = ".", x = genus)) %>% 
    dplyr::mutate(genus = gsub(pattern = "unclassified", replacement = "uncl.", x = genus)) %>% 
    dplyr::mutate(genus = gsub(pattern = "-", replacement = ".", x = genus)) %>% 
    dplyr::arrange(genus)

plot_data <- ps_16s_data %>% 
    microbiome::abundances(x = ., transform = "compositional") %>% 
    otu_table(taxa_are_rows = T) %>% t()%>% data.frame() 
plot_data$Diet <- ps_16s_data %>% sample_data() %>% data.frame() %>% .$Diet
colnames(plot_data) <- gsub(pattern = "unclassified", replacement = "uncl.", x = colnames(plot_data))

p_16s_healthy <- purrr::pmap(dplyr::select(da_16s_healthy_2, genus, score), function(genus, score) {
        ggplot(plot_data, aes_string("Diet", genus)) +
            geom_boxplot(aes(fill = Diet), outlier.shape = NA, fill = colv, width = 0.4) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            ggtitle(glue::glue("Score {score}")) +
        theme(axis.line = element_line(colour = "black"),
              legend.text = element_text(),
              legend.key = element_rect(fill = "transparent"),
              legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 0, size=12),
              #axis.text.x = element_blank(),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12, face = "italic"),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
})

p_16s_healthy_final <- wrap_plots(p_16s_healthy, ncol = 3) +
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Differential abundance analysis 16S rRNA healthy group') 

# 16S: DA DM II Group --------------------------------------------------------

ps_16s_data <- ps_16s_g %>% subset_samples(physeq = ., Disease == "DM II")

#
# ALDEx2
#

da_aldex <- aldex.clr(
    reads = ps_16s_data %>% otu_table(object = ., taxa_are_rows = T) %>% data.frame(),
    conds = ps_16s_data %>% sample_data(object = .) %>% data.frame() %>% .$Diet, 
    mc.samples = 512, 
    denom = "all",
    verbose = TRUE
)

# calculate expected values of the Welch's t-test and Wilcoxon rank test 
da_aldex_tt <- aldex.ttest(
    clr = da_aldex, 
    paired.test = FALSE, 
    verbose = TRUE)
# effect sizes
da_aldex_effect <- aldex.effect(da_aldex, CI = TRUE, verbose = FALSE)
# combine outputs 
aldex_16s_diab <- data.frame(da_aldex_tt, da_aldex_effect)

aldex_16s_diab %>% 
    rownames_to_column("genus") %>%
    dplyr::filter(we.eBH <= qcutoff)  %>% # t-test output
    dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>% 
    knitr::kable()

#
# ANCOMBC
#

# perform the analysis 
da_ancombc <- ancombc(
    phyloseq = ps_16s_data, 
    formula = "Diet", group = "Diet", tax_level = "Genus",
    p_adj_method = "BH", alpha = qcutoff, 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = TRUE, # recommended if the sample size is small 
    global = FALSE
)
# store the results in res 
da_ancombc_diab <- da_ancombc$res

da_ancombc_diab$diff_abn %>% 
    dplyr::filter(DietLFD == TRUE) %>% 
    knitr::kable()

#
# MaAsLin2
#

# similar settings as in Nearing et al. (2021)
da_maaslin <- Maaslin2( 
    input_data = ps_16s_data %>% otu_table(object = ., taxa_are_rows = T) %>% 
        data.frame() %>% t(), 
    input_metadata = ps_16s_data %>% sample_data() %>% data.frame(),
    output = "MaAsLin2_16s_diab",
    transform = "AST",
    fixed_effects = "Diet",
    reference = "Diet,HFD",  
    normalization = "TSS",
    standardize = FALSE, cores = 4, 
    plot_heatmap = T, heatmap_first_n = 10,
    min_prevalence = 0, # prev filterin already done
    correction = "BH"
)

da_maaslin$results %>% 
    dplyr::filter(qval < qcutoff)

#
# dacomp
#

df_dacomp <- ps_16s_data %>% 
    otu_table(object = ., taxa_are_rows = T) %>% 
    data.frame() %>% t()
result.selected.references <- dacomp.select_references(
    X = df_dacomp,
    minimal_TA = 100, 
    verbose = T)
length(result.selected.references$selected_references)
dacomp.plot_reference_scores(result.selected.references)
colnames(df_dacomp)[result.selected.references$selected_references]

#multiplicity correction levels for the BH and DS-FDR methods
q_BH <- q_DSFDR <- qcutoff

#Perform testing:
result.test <- dacomp.test(
    X = df_dacomp,
    y = ps_16s_data %>% 
        sample_data(object = .) %>% .$Diet,
    ind_reference_taxa = result.selected.references,
    test = DACOMP.TEST.NAME.WILCOXON,
    # nr_perm = 1000, # use default values here
    verbose = T,
    q = q_DSFDR,
    return_rarefied_values = TRUE)

rejected_BH <- which(p.adjust(result.test$p.values.test,method = 'BH') <= q_BH)
rejected_DSFDR <- result.test$dsfdr_rejected

da_dacomp <- data.frame(
    feature = df_dacomp %>% colnames(),
    pval = result.test$p.values.test,
    qval_DSFDR = result.test$p.values.test.adjusted,
    effect = result.test$effect_size_estimates
)

#
# Summarise methods
#
da_16s_diab <- dplyr::full_join(
    x = rownames_to_column(aldex_16s_diab, "genus") %>%
        dplyr::select(genus, aldex2 = we.eBH),
    y = rownames_to_column(da_ancombc$res$diff_abn, "genus") %>%
        dplyr::select(genus, ancombc = DietLFD),
    by = "genus") %>%
    dplyr::full_join( x = .,
                      y = dplyr::select(da_maaslin$results, genus = feature, maaslin2 = qval), 
                      by = "genus") %>%
    dplyr::full_join( x = .,
                      dplyr::select(da_dacomp, genus = feature, dacomp = qval_DSFDR), 
                      by = "genus") %>%
    dplyr::mutate(
        across(c(aldex2, maaslin2, dacomp), ~ .x <= qcutoff),
        score = rowSums(across(c(aldex2, ancombc, maaslin2, dacomp)), na.rm = T)
    ) %>% 
    dplyr::filter( !is.na(score) & score > 0)

da_16s_diab$aldex2[ is.na(da_16s_diab$aldex2) ] <- FALSE
da_16s_diab$ancombc[ is.na(da_16s_diab$ancombc) ] <- FALSE
da_16s_diab$maaslin2[ is.na(da_16s_diab$maaslin2) ] <- FALSE
da_16s_diab$dacomp[ is.na(da_16s_diab$dacomp) ] <- FALSE

# how many genera were identified by each method?
summarise(da_16s_diab, across(where(is.logical), sum)) %>%
    knitr::kable()

da_16s_diab_2 <- da_16s_diab %>% 
    dplyr::filter(score >= 1)

# ITS: DA Control Group --------------------------------------------------------

ps_ITS_data <- ps_its_g %>% subset_samples(physeq = ., Disease == "Healthy")

#
# ALDEx2
#

da_aldex <- aldex.clr(
    reads = ps_ITS_data %>% otu_table(object = ., taxa_are_rows = T) %>% data.frame(),
    conds = ps_ITS_data %>% sample_data(object = .) %>% .$Diet, 
    mc.samples = 512, 
    denom = "all",
    verbose = TRUE
)

# calculate expected values of the Welch's t-test and Wilcoxon rank test 
da_aldex_tt <- aldex.ttest(
    clr = da_aldex, 
    paired.test = FALSE, 
    verbose = TRUE)
# effect sizes
da_aldex_effect <- aldex.effect(da_aldex, CI = TRUE, verbose = FALSE)
# combine outputs 
aldex_ITS_healthy <- data.frame(da_aldex_tt, da_aldex_effect)

aldex_ITS_healthy %>% 
    rownames_to_column("genus") %>%
    dplyr::filter(we.eBH <= qcutoff)  %>% # t-test output
    dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>% 
    knitr::kable()

#
# ANCOMBC
#

# perform the analysis 
da_ancombc <- ancombc(
    phyloseq = ps_ITS_data, 
    formula = "Diet", group = "Diet", tax_level = "Genus",
    p_adj_method = "BH", alpha = qcutoff, 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = TRUE, # recommended if the sample size is small 
    global = FALSE
)
# store the results in res 
da_ancombc_healthy <- da_ancombc$res

da_ancombc_healthy$diff_abn %>% 
    dplyr::filter(DietLFD == TRUE) %>% 
    knitr::kable()

#
# MaAsLin2
#

# similar settings as in Nearing et al. (2021)
da_maaslin <- Maaslin2( 
    input_data = ps_ITS_data %>% otu_table(object = ., taxa_are_rows = T) %>% 
        data.frame() %>% t(), 
    input_metadata = ps_ITS_data %>% sample_data() %>% data.frame(),
    output = "MaAsLin2_ITS_healthy",
    transform = "AST",
    fixed_effects = "Diet",
    reference = "Diet,HFD",  
    normalization = "TSS",
    standardize = FALSE, cores = 4, 
    plot_heatmap = T, heatmap_first_n = 10,
    min_prevalence = 0, # prev filterin already done
    correction = "BH"
)

da_maaslin$results %>% 
    dplyr::filter(qval < qcutoff)

#
# dacomp
#

df_dacomp <- ps_ITS_data %>% 
    otu_table(object = ., taxa_are_rows = T) %>% 
    data.frame() %>% t()
result.selected.references <- dacomp.select_references(
    X = df_dacomp,
    minimal_TA = 100, 
    verbose = T)
length(result.selected.references$selected_references)
dacomp.plot_reference_scores(result.selected.references)
colnames(df_dacomp)[result.selected.references$selected_references]

#multiplicity correction levels for the BH and DS-FDR methods
q_BH <- q_DSFDR <- qcutoff

#Perform testing:
result.test <- dacomp.test(
    X = df_dacomp,
    y = ps_ITS_data %>% 
        sample_data(object = .) %>% .$Diet,
    ind_reference_taxa = result.selected.references,
    test = DACOMP.TEST.NAME.WILCOXON,
    # nr_perm = 1000, # use default values here
    verbose = T,
    q = q_DSFDR,
    return_rarefied_values = TRUE)

rejected_BH <- which(p.adjust(result.test$p.values.test,method = 'BH') <= q_BH)
rejected_DSFDR <- result.test$dsfdr_rejected

da_dacomp <- data.frame(
    feature = df_dacomp %>% colnames(),
    pval = result.test$p.values.test,
    qval_DSFDR = result.test$p.values.test.adjusted,
    effect = result.test$effect_size_estimates
)

#
# Summarise methods
#
da_ITS_healthy <- dplyr::full_join(
    x = rownames_to_column(aldex_ITS_healthy, "genus") %>%
        dplyr::select(genus, aldex2 = we.eBH),
    y = rownames_to_column(da_ancombc$res$diff_abn, "genus") %>%
        dplyr::select(genus, ancombc = DietLFD),
    by = "genus") %>%
    dplyr::full_join( x = .,
                      y = dplyr::select(da_maaslin$results, genus = feature, maaslin2 = qval), 
                      by = "genus") %>%
    dplyr::full_join( x = .,
                      dplyr::select(da_dacomp, genus = feature, dacomp = qval_DSFDR), 
                      by = "genus") %>%
    dplyr::mutate(
        across(c(aldex2, maaslin2, dacomp), ~ .x <= qcutoff),
        score = rowSums(across(c(aldex2, ancombc, maaslin2, dacomp)), na.rm = T)
    ) %>% 
    dplyr::filter( !is.na(score) & score > 0)
da_ITS_healthy$aldex2[ is.na(da_ITS_healthy$aldex2) ] <- FALSE
da_ITS_healthy$ancombc[ is.na(da_ITS_healthy$ancombc) ] <- FALSE
da_ITS_healthy$maaslin2[ is.na(da_ITS_healthy$maaslin2) ] <- FALSE
da_ITS_healthy$dacomp[ is.na(da_ITS_healthy$dacomp) ] <- FALSE

# how many genera were identified by each method?
summarise(da_ITS_healthy, across(where(is.logical), sum)) %>%
    knitr::kable()

da_ITS_healthy_2 <- da_ITS_healthy %>% 
    dplyr::filter(score >= 2) %>% 
    dplyr::mutate(genus = gsub(pattern = "unclassified", replacement = "uncl.", x = genus))

# ITS: DA DM II Group --------------------------------------------------------

ps_ITS_data <- ps_its_g %>% subset_samples(physeq = ., Disease == "DM II")

#
# ALDEx2
#

da_aldex <- aldex.clr(
    reads = ps_ITS_data %>% otu_table(object = ., taxa_are_rows = T) %>% data.frame(),
    conds = ps_ITS_data %>% sample_data(object = .) %>% data.frame() %>% .$Diet, 
    mc.samples = 512, 
    denom = "all",
    verbose = TRUE
)

# calculate expected values of the Welch's t-test and Wilcoxon rank test 
da_aldex_tt <- aldex.ttest(
    clr = da_aldex, 
    paired.test = FALSE, 
    verbose = TRUE)
# effect sizes
da_aldex_effect <- aldex.effect(da_aldex, CI = TRUE, verbose = FALSE)
# combine outputs 
aldex_ITS_diab <- data.frame(da_aldex_tt, da_aldex_effect)

aldex_ITS_diab %>% 
    rownames_to_column("genus") %>%
    dplyr::filter(we.eBH <= qcutoff)  %>% # t-test output
    dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>% 
    knitr::kable()

#
# ANCOMBC
#

# perform the analysis 
da_ancombc <- ancombc(
    phyloseq = ps_ITS_data, 
    formula = "Diet", group = "Diet", tax_level = "Genus",
    p_adj_method = "BH", alpha = qcutoff, 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = TRUE, # recommended if the sample size is small 
    global = FALSE
)
# store the results in res 
da_ancombc_diab <- da_ancombc$res

da_ancombc_diab$diff_abn %>% 
    dplyr::filter(DietLFD == TRUE) %>% 
    knitr::kable()

#
# MaAsLin2
#

# similar settings as in Nearing et al. (2021)
da_maaslin <- Maaslin2( 
    input_data = ps_ITS_data %>% otu_table(object = ., taxa_are_rows = T) %>% 
        data.frame() %>% t(), 
    input_metadata = ps_ITS_data %>% sample_data() %>% data.frame(),
    output = "MaAsLin2_ITS_diab",
    transform = "AST",
    fixed_effects = "Diet",
    reference = "Diet,HFD",  
    normalization = "TSS",
    standardize = FALSE, cores = 4, 
    plot_heatmap = T, heatmap_first_n = 10,
    min_prevalence = 0, # prev filterin already done
    correction = "BH"
)

da_maaslin$results %>% 
    dplyr::filter(qval < qcutoff)

#
# dacomp
#

df_dacomp <- ps_ITS_data %>% 
    otu_table(object = ., taxa_are_rows = T) %>% 
    data.frame() %>% t()
result.selected.references <- dacomp.select_references(
    X = df_dacomp,
    minimal_TA = 100, 
    verbose = T)
length(result.selected.references$selected_references)
dacomp.plot_reference_scores(result.selected.references)
colnames(df_dacomp)[result.selected.references$selected_references]

#multiplicity correction levels for the BH and DS-FDR methods
q_BH <- q_DSFDR <- qcutoff

#Perform testing:
result.test <- dacomp.test(
    X = df_dacomp,
    y = ps_ITS_data %>% 
        sample_data(object = .) %>% .$Diet,
    ind_reference_taxa = result.selected.references,
    test = DACOMP.TEST.NAME.WILCOXON,
    # nr_perm = 1000, # use default values here
    verbose = T,
    q = q_DSFDR,
    return_rarefied_values = TRUE)

rejected_BH <- which(p.adjust(result.test$p.values.test,method = 'BH') <= q_BH)
rejected_DSFDR <- result.test$dsfdr_rejected

da_dacomp <- data.frame(
    feature = df_dacomp %>% colnames(),
    pval = result.test$p.values.test,
    qval_DSFDR = result.test$p.values.test.adjusted,
    effect = result.test$effect_size_estimates
)


# Summarise methods
da_ITS_diab <- dplyr::full_join(
    x = rownames_to_column(aldex_ITS_healthy, "genus") %>%
        dplyr::select(genus, aldex2 = we.eBH),
    y = rownames_to_column(da_ancombc$res$diff_abn, "genus") %>%
        dplyr::select(genus, ancombc = DietLFD),
    by = "genus") %>%
    dplyr::full_join( x = .,
                      y = dplyr::select(da_maaslin$results, genus = feature, maaslin2 = qval), 
                      by = "genus") %>%
    dplyr::full_join( x = .,
                      dplyr::select(da_dacomp, genus = feature, dacomp = qval_DSFDR), 
                      by = "genus") %>%
    dplyr::mutate(
        across(c(aldex2, maaslin2, dacomp), ~ .x <= qcutoff),
        score = rowSums(across(c(aldex2, ancombc, maaslin2, dacomp)), na.rm = T)
    ) %>% 
    dplyr::filter( !is.na(score) & score > 0)
da_ITS_diab$aldex2[ is.na(da_ITS_diab$aldex2) ] <- FALSE
da_ITS_diab$ancombc[ is.na(da_ITS_diab$ancombc) ] <- FALSE
da_ITS_diab$maaslin2[ is.na(da_ITS_diab$maaslin2) ] <- FALSE
da_ITS_diab$dacomp[ is.na(da_ITS_diab$dacomp) ] <- FALSE

# how many genera were identified by each method?
summarise(da_ITS_diab, across(where(is.logical), sum)) %>%
    knitr::kable()

da_ITS_diab_2 <- da_ITS_diab %>% 
    dplyr::filter(score >= 1)

# Metabolome: DA Healthy Group --------------------------------------------------------

ps_met_data <- ps_met_g %>% subset_samples(physeq = ., Disease == "Healthy")

#
# ALDEx2
#

da_aldex <- aldex.clr(
    reads = ps_met_data %>% otu_table(object = ., taxa_are_rows = T) %>% t() %>% data.frame() %>% ceiling(),
    conds = ps_met_data %>% sample_data(object = .) %>% .$Diet, 
    mc.samples = 512, 
    denom = "all",
    verbose = TRUE
)

# calculate expected values of the Welch's t-test and Wilcoxon rank test 
da_aldex_tt <- aldex.ttest(
    clr = da_aldex, 
    paired.test = FALSE, 
    verbose = TRUE)
# effect sizes
da_aldex_effect <- aldex.effect(da_aldex, CI = TRUE, verbose = FALSE)
# combine outputs 
aldex_met_healthy <- data.frame(da_aldex_tt, da_aldex_effect)

aldex_met_healthy %>% 
    rownames_to_column("metabolite") %>%
    dplyr::filter(we.eBH <= qcutoff)  %>% # t-test output
    dplyr::select(metabolite, we.eBH, wi.eBH, effect, overlap) %>% 
    knitr::kable()

#
# ANCOMBC
#

# perform the analysis 
da_ancombc <- ancombc(
    phyloseq = ps_met_data, 
    formula = "Diet", group = "Diet", 
    p_adj_method = "BH", alpha = qcutoff, 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = TRUE, # recommended if the sample size is small 
    global = FALSE
)
# store the results in res 
da_ancombc_healthy <- da_ancombc$res

da_ancombc_healthy$diff_abn %>% 
    dplyr::filter(DietLFD == TRUE) %>% 
    knitr::kable()

#
# MaAsLin2
#

# similar settings as in Nearing et al. (2021)
da_maaslin <- Maaslin2( 
    input_data = ps_met_data %>% otu_table(object = ., taxa_are_rows = T) %>% 
        data.frame(), 
    input_metadata = ps_met_data %>% sample_data() %>% data.frame(),
    output = "MaAsLin2_Metabolome_healthy",
    transform = "AST",
    fixed_effects = "Diet",
    reference = "Diet,HFD",  
    normalization = "TSS",
    standardize = FALSE, cores = 4, 
    plot_heatmap = T, heatmap_first_n = 10,
    min_prevalence = 0, # prev filterin already done
    correction = "BH"
)

da_maaslin$results %>% 
    dplyr::filter(qval < qcutoff)

#
# dacomp --> not possible because it runs on counts and numbers are too low
#

#
# Summarise methods
#
da_met_healthy <- dplyr::full_join(
    rownames_to_column(aldex_met_healthy, "metabolite") %>%
        dplyr::select(metabolite, aldex2 = we.eBH),
    rownames_to_column(da_ancombc$res$diff_abn, "metabolite") %>%
        dplyr::select(metabolite, ancombc = DietLFD),
    by = "metabolite") %>%
    dplyr::full_join(
        dplyr::select(da_maaslin$results, metabolite = feature, maaslin2 = qval), 
        by = "metabolite") %>%
    dplyr::mutate(
        across(c(aldex2, maaslin2), ~ .x <= qcutoff),
        score = rowSums(across(c(aldex2, ancombc, maaslin2)), na.rm = T)
    ) %>% 
    dplyr::filter( !is.na(score) & score > 0)

da_met_healthy$aldex2[ is.na(da_met_healthy$aldex2) ] <- FALSE
da_met_healthy$ancombc[ is.na(da_met_healthy$ancombc) ] <- FALSE
da_met_healthy$maaslin2[ is.na(da_met_healthy$maaslin2) ] <- FALSE

# how many genera were identified by each method?
summarise(da_met_healthy, across(where(is.logical), sum)) %>%
    knitr::kable()

da_met_healthy_2 <- da_met_healthy %>% 
    dplyr::filter(score >= 2) %>% 
    dplyr::arrange(metabolite)
da_met_healthy_2$metabolite[da_met_healthy_2$metabolite == "Proprionate"] <- "Propionate"

plot_data <- ps_met_data %>% 
    # microbiome::abundances(x = ., transform = "compositional") %>% 
    otu_table(taxa_are_rows = T) %>% 
    # t() %>% 
    data.frame() %>% 
    rename(Propionate = Proprionate) 
plot_data$Diet <- ps_met_data %>% sample_data() %>% data.frame() %>% .$Diet

p_met_healthy <- purrr::pmap(dplyr::select(da_met_healthy_2, metabolite, score), function(metabolite, score) {
    ggplot(plot_data, aes_string("Diet", metabolite)) +
        geom_boxplot(aes(fill = Diet), outlier.shape = NA, fill = colv, width = 0.4) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        ggtitle(glue::glue("Score {score}")) +
        theme(axis.line = element_line(colour = "black"),
              legend.text = element_text(),
              legend.key = element_rect(fill = "transparent"),
              legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 0, size=12),
              #axis.text.x = element_blank(),
              axis.text.y = element_text(size=12),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
})

p_met_healthy_final <- wrap_plots(p_met_healthy) +
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Differential abundance analysis Metabolites healthy group') 

# Metabolome: DA DM II Group --------------------------------------------------------

ps_met_data <- ps_met_g %>% subset_samples(physeq = ., Disease == "DM II")

#
# ALDEx2
#

da_aldex <- aldex.clr(
    reads = ps_met_data %>% otu_table(object = ., taxa_are_rows = T) %>% t() %>% data.frame() %>% ceiling(),
    conds = ps_met_data %>% sample_data(object = .) %>% .$Diet, 
    mc.samples = 512, 
    denom = "all",
    verbose = TRUE
)

# calculate expected values of the Welch's t-test and Wilcoxon rank test 
da_aldex_tt <- aldex.ttest(
    clr = da_aldex, 
    paired.test = FALSE, 
    verbose = TRUE)
# effect sizes
da_aldex_effect <- aldex.effect(da_aldex, CI = TRUE, verbose = FALSE)
# combine outputs 
aldex_met_disease <- data.frame(da_aldex_tt, da_aldex_effect)

aldex_met_disease %>% 
    rownames_to_column("metabolite") %>%
    dplyr::filter(we.eBH <= qcutoff)  %>% # t-test output
    dplyr::select(metabolite, we.eBH, wi.eBH, effect, overlap) %>% 
    knitr::kable()

#
# ANCOMBC
#

# perform the analysis 
da_ancombc <- ancombc(
    phyloseq = ps_met_data, 
    formula = "Diet", group = "Diet", 
    p_adj_method = "BH", alpha = qcutoff, 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = TRUE, # recommended if the sample size is small 
    global = FALSE
)
# store the results in res 
da_ancombc_disease <- da_ancombc$res

da_ancombc_disease$diff_abn %>% 
    dplyr::filter(DietLFD == TRUE) %>% 
    knitr::kable()

#
# MaAsLin2
#

# similar settings as in Nearing et al. (2021)
da_maaslin <- Maaslin2( 
    input_data = ps_met_data %>% otu_table(object = ., taxa_are_rows = T) %>% 
        data.frame(), 
    input_metadata = ps_met_data %>% sample_data() %>% data.frame(),
    output = "MaAsLin2_Metabolome_disease",
    transform = "AST",
    fixed_effects = "Diet",
    reference = "Diet,HFD",  
    normalization = "TSS",
    standardize = FALSE, cores = 4, 
    plot_heatmap = T, heatmap_first_n = 10,
    min_prevalence = 0, # prev filterin already done
    correction = "BH"
)

da_maaslin$results %>% 
    dplyr::filter(qval < qcutoff)

# Summarise methods
da_met_disease <- dplyr::full_join(
    rownames_to_column(aldex_met_disease, "metabolite") %>%
        dplyr::select(metabolite, aldex2 = we.eBH),
    rownames_to_column(da_ancombc$res$diff_abn, "metabolite") %>%
        dplyr::select(metabolite, ancombc = DietLFD),
    by = "metabolite") %>%
    dplyr::full_join(
        dplyr::select(da_maaslin$results, metabolite = feature, maaslin2 = qval), 
        by = "metabolite") %>%
    dplyr::mutate(
        across(c(aldex2, maaslin2), ~ .x <= qcutoff),
        score = rowSums(across(c(aldex2, ancombc, maaslin2)), na.rm = T)
    ) %>% 
    dplyr::filter( !is.na(score) & score > 0)


# 16S, ITS, Metabolome: Summary -----------------------------------------------------------------

# 16S
summarise(da_16s_healthy, across(where(is.logical), sum)) %>%
    knitr::kable()
da_16s_healthy %>% 
    dplyr::filter(score >= 2)

summarise(da_16s_diab, across(where(is.logical), sum)) %>%
    knitr::kable()
da_16s_diab %>% 
    dplyr::filter(score >= 2)

# ITS
summarise(da_ITS_healthy, across(where(is.logical), sum)) %>%
    knitr::kable()
da_ITS_healthy %>% 
    dplyr::filter(score >= 2)

summarise(da_ITS_diab, across(where(is.logical), sum)) %>%
    knitr::kable()
da_ITS_diab %>% 
    dplyr::filter(score >= 2)

# Metabolome
summarise(da_met_healthy, across(where(is.logical), sum)) %>%
    knitr::kable()
da_met_healthy %>% 
    dplyr::filter(score >= 2)

summarise(da_met_disease, across(where(is.logical), sum)) %>%
    knitr::kable()
da_met_disease %>% 
    dplyr::filter(score >= 2)

# Pretty plotting ---------------------------------------------------------

p_16s_healthy_final <- wrap_plots(p_16s_healthy, ncol = 3) +
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'A')
ggsave(filename = "plots_final/Suppl_Fig_4.pdf", width = 9, height = 9, plot = p_16s_healthy_final)

p_met_healthy_final <- wrap_plots(p_met_healthy) +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = 'A')
ggsave(filename = "plots_final/Suppl_Fig_5.pdf", width = 6, height = 6, plot = p_met_healthy_final)

my_list <- list()
my_list[["Bacterial control fdr < 0.1"]] <- da_16s_healthy
my_list[["Bacterial DM II fdr < 0.1"]] <- da_16s_diab
my_list[["Fungal DM II fdr < 0.1"]] <- da_ITS_healthy
my_list[["Fungal healthy fdr < 0.1"]] <- da_ITS_diab
my_list[["Metabolites control fdr < 0.1"]] <- da_met_healthy
my_list[["Metabolites DM II fdr < 0.1"]] <- da_met_disease

WriteXLS::WriteXLS(x = my_list, ExcelFileName = "Suppl_Table_DA.xlsx", 
                   row.names = F, 
                   AdjWidth = T, 
                   FreezeRow = 1, 
                   BoldHeaderRow = TRUE)
