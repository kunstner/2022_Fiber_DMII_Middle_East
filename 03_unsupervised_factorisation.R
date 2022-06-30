
# MOFA --------------------------------------------------------------------

# https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/microbiome_vignette.html
# https://github.com/bwhaak/MOFA_microbiome

# https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)

library(ggpubr)
library(patchwork)

library(MOFA2)
library(compositions)

# Data --------------------------------------------------------------------

ps_16s <- readRDS(file = "data/processed_ps_16s.RDS")
ps_its <- readRDS(file = "data/processed_ps_its.RDS")
metabo <- readRDS(file = "data/processed_metabo.RDS")
cova   <- readRDS(file = "data/processed_cova.RDS")

seed <- 42
colv <- c("steelblue", "orange")
colv_groups <- c("forestgreen", "mediumpurple4" )

# Format data -------------------------------------------------------------

cova <- cova %>% 
    rename(Diet = Fiber.Content.Diet) %>% 
    rename(Weight = Weight.KG.) %>% 
    rename(BMI = BMI..Kg.m) 

cova$group <- "DM II"
cova$group[grep(pattern = "control", x = cova$Status) ] <- "Control"

ps_16s <- ps_16s %>% tax_glom(physeq = ., taxrank = "Genus")
ps_its <- ps_its %>% tax_glom(physeq = ., taxrank = "Genus")

bact <- otu_table(object = ps_16s) %>% t() %>% data.frame()
#colnames(bact) <- data.frame( tax_table( object = ps_16s ) )[, "Genus"]

fung <- otu_table(object = ps_its) %>% t() %>% data.frame()
#colnames(fung) <- data.frame( tax_table( object = ps_its ) )[, "Genus"]

rownames(bact) == rownames(fung) 

metabo <- metabo %>% 
    column_to_rownames("Sample") 

# remove taxa with low variability
bact <- bact[ , apply(bact, 2, var) > 10 ]
fung <- fung[ , apply(fung, 2, var) > 10 ]
metabo <- metabo[ , apply(metabo, 2, var) > 5 ]

ps_16s <- prune_taxa(taxa = colnames(bact), x = ps_16s) 
ps_its <- prune_taxa(taxa = colnames(fung), x = ps_its) 

# Normalise with CLR

norm_b.clr <- microbiome::transform(x = ps_16s, transform = "clr") %>% 
    otu_table() %>% t() %>% data.frame()
colnames(norm_b.clr) <- data.frame( tax_table( object = ps_16s ) )[, "Genus"]
norm_b.clr <- norm_b.clr %>% 
    rownames_to_column(var = "sample") %>%
    reshape::melt(id.var = "sample", variable.name = "feature") %>%
    dplyr::mutate(view = "Bacteria") %>% 
    rename(feature = variable) %>% 
    dplyr::mutate(feature = gsub(pattern = "g__", replacement = "", x = feature)) %>% 
    dplyr::mutate(feature = gsub(pattern = "unclassified", replacement = "uncl", x = feature)) 

norm_f.clr <- microbiome::transform(x = ps_its, transform = "clr") %>% 
    otu_table() %>% t() %>% data.frame()
colnames(norm_f.clr) <- data.frame( tax_table( object = ps_its ) )[, "Genus"]
norm_f.clr <- norm_f.clr %>% 
    rownames_to_column(var = "sample") %>%
    reshape::melt(id.var = "sample", variable.name = "feature") %>%
    dplyr::mutate(view = "Fungi") %>% 
    rename(feature = variable) %>% 
    dplyr::mutate(feature = gsub(pattern = "unclassified", replacement = "uncl", x = feature)) 

norm_m.clr <- as.data.frame( compositions::clr(metabo+0.1)) %>%
    rownames_to_column(var = "sample") %>%
    reshape::melt(id.var = "sample", variable.name = "feature") %>%
    dplyr::mutate(view = "Metabolome") %>% 
    rename(feature = variable) 

saveRDS(object = norm_b.clr, file = "data/norm_b.clr.RDS")
saveRDS(object = norm_f.clr, file = "data/norm_f.clr.RDS")
saveRDS(object = norm_m.clr, file = "data/norm_m.clr.RDS")

norm_b.clr$group <- "DM II"
norm_b.clr$group[grep(pattern = "C", x = norm_b.clr$sample) ] <- "Control"

norm_f.clr$group <- "DM II"
norm_f.clr$group[grep(pattern = "C", x = norm_f.clr$sample) ] <- "Control"

norm_m.clr$group <- "DM II"
norm_m.clr$group[grep(pattern = "C", x = norm_m.clr$sample) ] <- "Control"

data <- rbind(norm_b.clr, norm_f.clr, norm_m.clr)
head( data )

unique(data$view)

ggdensity(data, x="value", fill="gray70") +
    facet_wrap(~view + group, nrow=3, scales="free") 

# Prepare MOFA object and training ----------------------------------------

MOFAobject <- create_mofa_from_df( df = data )

p_mofa_over <- plot_data_overview(object = MOFAobject, colors = c(pals::tableau20(20)[c(2,4,6)]) )
p_mofa_over

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- seed
train_opts$maxiter <- 10000 # default 1,000

# Train model
MOFAobject <- prepare_mofa(object = MOFAobject, 
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts,
)

# Run MOFA and add Metadata -----------------------------------------------

MOFAobject <- run_mofa(object = MOFAobject, 
                       outfile="MOFA2/MOFA2_trained_diabetes.hdf5", 
                       use_basilisk = TRUE)
# Iteration 3286: time=0.01, ELBO=-14674.54, deltaELBO=-0.000 (0.00000036%), Factors=10
# Converged!
MOFAobject

# saveRDS(object = MOFAobject, file = "MOFAobject.RDS")
MOFAobject <- readRDS(file = "data/MOFAobject.RDS")

p_mofa_over <- plot_data_overview(object = MOFAobject, colors = c(pals::tableau20(20)[c(2,4,6)]) )
p_mofa_over

# add metadata
samples_metadata(MOFAobject) <- cova %>% 
    rename(sample = Sample) # %>% rename(group = Status)

plot_factor_cor(object = MOFAobject)

# variance decomposition
calculate_variance_explained(object = MOFAobject, groups = 'all')

p_mofa_variance1 <- plot_variance_explained(MOFAobject, plot_total = T)[[2]] +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12))
p_mofa_variance1

p_mofa_variance2 <- plot_variance_explained(MOFAobject, max_r2=50) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12))
p_mofa_variance2

r2 <- MOFAobject@cache$variance_explained$r2_per_factor[[1]]
r2.dt <- r2 %>%
    data.frame() %>% 
    dplyr::mutate(factor = as.factor(1:MOFAobject@dimensions$K)) %>% 
    reshape::melt(id.vars=c("factor"), variable.name="view", value.name = "r2") %>%
    dplyr::rename(view = variable) 
r2.dt$cum_r2[r2.dt$view == "Bacteria"] <- cumsum(x = r2.dt$value[r2.dt$view == "Bacteria"])
r2.dt$cum_r2[r2.dt$view == "Fungi"] <- cumsum(x = r2.dt$value[r2.dt$view == "Fungi"])
r2.dt$cum_r2[r2.dt$view == "Metabolome"] <- cumsum(x = r2.dt$value[r2.dt$view == "Metabolome"])

ggline(r2.dt, x="factor", y="cum_r2", color="view") +
    labs(x="Factor number", y="Cumulative variance explained (%)") +
    theme(
        legend.title = element_blank(), 
        legend.position = "top",
        axis.text = element_text(size=rel(0.8))
    )

# Association analysis
p_cov_all <- correlate_factors_with_covariates(object = MOFAobject, transpose = T, 
                                  covariates = c("Diet", "Gender", "Age", "BMI"), 
                                  plot="log_pval", alpha = 0.01 ) 

p_cov_ctr <- correlate_factors_with_covariates(object = MOFAobject, groups = "Control", transpose = T, 
                                  covariates = c("Diet", "Gender", "Age", "BMI"), 
                                       plot="log_pval", alpha = 0.01 ) 

p_cov_dia <- correlate_factors_with_covariates(object = MOFAobject, groups = "DM II", transpose = T,
                                  covariates = c("Diet", "Gender", "Age", "BMI"), 
                                  # return_data = TRUE,
                                  plot="log_pval", alpha = 0.01 ) 

cowplot::plot_grid(p_cov_all$gtable, p_cov_ctr$gtable)

# Diet ----------------------------------------------------------------

p_diet <- plot_factors(MOFAobject, 
                  factors = c(1), 
                  color_by = "Diet", shape_by = "group",
                  dot_size = 4) + 
    scale_fill_manual(values=colv) + 
#    geom_density_2d(aes_string(color="color_by")) +
    stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.05, level = 0.95) + 
    scale_color_manual(values=colv) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          #axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_diet

p_factor_vio <- plot_factor(MOFAobject, 
            factor = c(1), 
            color_by = "Diet", shape_by = "group",
            dot_size = 2,
            dodge = TRUE,
            stroke = 0.4, 
            add_violin = T, violin_alpha = 0.25,
            add_boxplot = F, boxplot_alpha = 0.25) +
    scale_fill_manual(values=colv) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=0),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=12),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid"))
p_factor_vio

# Characterization of Factor 1 --------------------------------------------

# Bacteria
p_f1_top_weights_bact <- plot_top_weights(object = MOFAobject, 
             factors = c(1), 
             view = "Bacteria", scale = TRUE,
             nfeatures = 6)

p_f1_weights_bact <- plot_weights(object = MOFAobject, 
             factors = c(1), view = "Bacteria", 
             nfeatures = 6)

p_f1_scatter_bact <- plot_data_scatter(object = MOFAobject, 
                  view = "Bacteria",
                  factor = 1, 
                  shape_by = "group",
                  features = 6, alpha = 1, 
                  add_lm = T, lm_per_group = T, dot_size = 2,
                  sign = "all",
                  color_by = "Diet" ) + 
    ylab("CLR transformed bacterial abundances") +
    scale_fill_manual(values=colv) +
    scale_color_manual(values=colv_groups) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=12)) +
    ylim(-10, 10) 
p_f1_scatter_bact

p_heat_bact <- plot_data_heatmap(MOFAobject, 
                  factor = 1, 
                  view = "Bacteria", 
                  features = 6,
                  denoise = TRUE,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  annotation_samples = c("Diet", "group"),  
                  annotation_colors = list(
                      Diet = c("High fiber diet"=colv[1], "Low fiber diet"=colv[2]),
                      group = c("Control"=colv_groups[1], "DM II"=colv_groups[2])), 
                  annotation_legend = F,
                  scale = "row", silent = TRUE
)

p_heat_bact <- ggplotify::as.ggplot(p_heat_bact)
p_f1_top_weights_bact
p_f1_weights_bact
p_f1_scatter_bact

layout <- "
AABBBB
CCCCCC
CCCCCC
"
p_f1_top_weights_bact + p_heat_bact +
    p_f1_scatter_bact +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Top Bacterial Contributions to Factor1') 

# Fungi
p_f1_top_weights_fung <- plot_top_weights(object = MOFAobject, 
                                          factor = 1, 
                                          view = "Fungi", 
                                          nfeatures = 3)

p_f1_weights_fung <- plot_weights(object = MOFAobject, 
                                  factor = 1, view = "Fungi", 
                                  nfeatures = 3)

p_f1_scatter_fung <- plot_data_scatter(MOFAobject, 
                                       view = "Fungi",
                                       factor = 1, 
                                       shape_by = "group",
                                       features = 3, alpha = 1, 
                                       add_lm = T, lm_per_group = T, dot_size = 2,
                                       sign = "all",
                                       color_by = "Diet" ) + 
    ylab("CLR transformed fungal abundances") +
    scale_fill_manual(values=colv) +
    scale_color_manual(values=colv_groups) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=12)) +
    ylim(-10, 10)

p_heat_fung <- plot_data_heatmap(MOFAobject, 
                                 factor = 1, 
                                 view = "Fungi", 
                                 features = 3,
                                 denoise = TRUE,
                                 cluster_rows = T, cluster_cols = F,
                                 show_colnames = F, show_rownames = T,
                                 annotation_samples = c("Diet", "group"),  
                                 annotation_colors = list(
                                     Diet = c("High fiber diet"=colv[1], "Low fiber diet"=colv[2]),
                                     group = c("Control"=colv_groups[1], "DM II"=colv_groups[2])), 
                                 annotation_legend = F,
                                 scale = "row", silent = TRUE
)

p_heat_fung <- ggplotify::as.ggplot(p_heat_fung)
p_f1_top_weights_fung
p_f1_weights_fung
p_f1_scatter_fung

layout <- "
AABBBB
CCCCCC
"
p_f1_top_weights_fung + p_heat_fung +
    p_f1_scatter_fung +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Top Fungal Contributions to Factor1') 

# Metabolome
p_f1_top_weights_meta <- plot_top_weights(object = MOFAobject, 
                                          factor = 1, 
                                          view = "Metabolome", 
                                          nfeatures = 4)

p_f1_weights_meta <- plot_weights(object = MOFAobject, 
                                  factor = 1, view = "Metabolome", 
                                  nfeatures = 4)

p_f1_scatter_meta <- plot_data_scatter(MOFAobject, 
                                       view = "Metabolome",
                                       factor = 1, 
                                       shape_by = "group",
                                       features = 4, alpha = 1, 
                                       add_lm = T, lm_per_group = T, dot_size = 2,
                                       sign = "all",
                                       color_by = "Diet" ) + 
    ylab("CLR transformed metabolite abundances") +
    scale_fill_manual(values=colv) +
    scale_color_manual(values=colv_groups) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=12)) +
    ylim(-10, 10)

p_heat_meta <- plot_data_heatmap(MOFAobject, 
                                 factor = 1, 
                                 view = "Metabolome", 
                                 features = 4,
                                 denoise = TRUE,
                                 cluster_rows = T, cluster_cols = F,
                                 show_colnames = F, show_rownames = T,
                                 annotation_samples = c("Diet", "group"),  
                                 annotation_colors = list(
                                     Diet = c("High fiber diet"=colv[1], "Low fiber diet"=colv[2]),
                                     group = c("Control"=colv_groups[1], "DM II"=colv_groups[2])), 
                                 annotation_legend = F, 
                                 scale = "row", silent = TRUE
)

p_heat_meta <- ggplotify::as.ggplot(p_heat_meta)
p_f1_top_weights_meta
p_f1_weights_meta
p_f1_scatter_meta

layout <- "
AABBBB
CCCCCC
"
p_f1_top_weights_meta + p_heat_meta +
    p_f1_scatter_meta +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Top Metabolome Contributions to Factor1') 

# Pretty plotting ---------------------------------------------------------

# Fig 3
layout <- "
ABBCC
ABBCC
"
p_mofa_over +
    p_factor_vio + p_mofa_variance2 +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A')
ggsave(filename = "plots_final/Figure_3.pdf", height = 6, width = 12)

# Suppl Fig 9
cowplot::plot_grid(p_cov_all$gtable, p_cov_ctr$gtable, nrow = 1, labels = c('A', 'B'))
ggsave(filename = "plots_final/Suppl_Fig_9.pdf", height = 4.5, width = 9.0)

# Suppl Fig 10
layout <- "
AABBBBBB
"
p_f1_top_weights_bact + 
    p_f1_scatter_bact +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Suppl_Fig_10.pdf", height = 6, width = 12)

# Suppl Fig 11
layout <- "
AABBBBBB
"
p_f1_top_weights_fung + 
    p_f1_scatter_fung +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Suppl_Fig_11.pdf", height = 4, width = 12)

# Suppl Fig 12
layout <- "
AABBBB
XXBBBB
"
p_f1_top_weights_meta + 
    p_f1_scatter_meta +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') 
ggsave(filename = "plots_final/Suppl_Fig_12.pdf", height = 6, width = 9)
