
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)

library(gbm)

# Data --------------------------------------------------------------------

seed <- 42
set.seed(seed)

norm_b.clr <- readRDS(file = "data/norm_b.clr.RDS") %>% 
    dplyr::select(-view) %>% 
    dplyr::mutate(feature = gsub(pattern = "\\[", replacement = "x", x = feature) ) %>% 
    dplyr::mutate(feature = gsub(pattern = "\\]", replacement = "x", x = feature) ) %>% 
    pivot_wider(data = ., names_from = feature, values_from = value) %>% 
    column_to_rownames('sample')

norm_f.clr <- readRDS(file = "data/norm_f.clr.RDS") %>% 
    dplyr::select(-view) %>% 
    pivot_wider(data = ., names_from = feature, values_from = value) %>% 
    column_to_rownames('sample')

norm_m.clr <- readRDS(file = "data/norm_m.clr.RDS") %>% 
    dplyr::select(-view) %>% 
    pivot_wider(data = ., names_from = feature, values_from = value) %>% 
    column_to_rownames('sample')

cova   <- readRDS(file = "data/processed_cova.RDS")

rownames(norm_b.clr) == cova$Sample
rownames(norm_f.clr) == cova$Sample
rownames(norm_m.clr) == cova$Sample

colnames(norm_b.clr) <- paste0("B_", colnames(norm_b.clr))
colnames(norm_f.clr) <- paste0("F_", colnames(norm_f.clr))
colnames(norm_m.clr) <- paste0("M_", colnames(norm_m.clr))

data <- cbind(norm_b.clr, norm_f.clr, norm_m.clr)

cova <- cova %>% 
    dplyr::filter(Status == "control") %>% 
    rename(Diet = Fiber.Content.Diet) %>% 
    rename(Weight = Weight.KG.) %>% 
    rename(BMI = BMI..Kg.m) %>% 
    dplyr::mutate(Diet = gsub(pattern = "High fiber diet", replacement = "HFD", x = Diet)) %>% 
    dplyr::mutate(Diet = gsub(pattern = "Low fiber diet", replacement = "LFD", x = Diet)) 

data <- data[ rownames(data) %in% cova$Sample, ]

data$Group <- cova$Diet

# 1,000-fold Gradient boosting --------------------------------------------

folds <- 1000
max_trees <- 20000

# re-format data for Bernoulli coding
data_folds <- data
data_folds$Group[data_folds$Group == "HFD"] <- 0
data_folds$Group[data_folds$Group == "LFD"] <- 1
data_folds$Group <- data_folds$Group %>% as.numeric

train.index <- caret::createDataPartition(data_folds$Group, p = .80, list = FALSE, times = folds)
cm_fold <- data.frame(i = 0, 
                      Accuracy = 0, 
                      Kappa = 0, 
                      AccuracyLower = 0, 
                      AccuracyUpper = 0, AccuracyNull = 0, 
                      AccuracyPValue = 0, 
                      McnemarPValue = 0,
                      BestTree = 0)

gmp_fold_rel <- data.frame(row.names = sort( colnames(data_folds %>% dplyr::select(-Group)) ))

for( i in 1:folds ) {
    print(paste0('Working on fold ', i))
    
    train_folds <- data_folds[ train.index[,i],]
    test_folds  <- data_folds[-train.index[,i],]
    # train_folds$Group <- factor(train_folds$Group)
    # test_folds$Group  <- factor(test_folds$Group)
    
    # run gradient boosting
    mod_gbm_fold <- gbm::gbm( Group ~.,
                        data = train_folds,
                        distribution = "bernoulli",
                        cv.folds = 5, 
                        class.stratify.cv = TRUE,
                        n.cores = 10,
                        shrinkage = .001, 
                        n.minobsinnode = 2, verbose = F,
                        n.trees = max_trees)
    # prediction step
    pred_fold <- gbm::predict.gbm(object = mod_gbm_fold,
                             newdata = test_folds,
                             n.trees = max_trees,
                             type = "response")
    pred_fold <- ifelse( test = pred_fold > 0.5, yes = 1, no = 0 ) # format predictions
    
    # confusion matrix
    cm <- caret::confusionMatrix(factor(test_folds$Group, levels = c(0,1)), 
                          factor(pred_fold, levels = c(0,1)) )
    cm_fold <- rbind(cm_fold, c(i, cm$overall) )
    
    best_tree <- gbm::gbm.perf(object = mod_gbm_fold, plot.it = F)
    cm_fold$BestTree[i] <-  best_tree
    
    gbm_imp <- gbm::summary.gbm(object = mod_gbm_fold, n.trees = best_tree, 
                                plotit = F, cBars = 5, las = 2)
    gbm_imp <- gbm_imp[sort( colnames(data_folds %>% dplyr::select(-Group)) ), ]
    gmp_fold_rel <- cbind(gmp_fold_rel, gbm_imp$rel.inf)
}

cm_fold <- cm_fold %>% 
    dplyr::filter(i != 0)

saveRDS(object = cm_fold, file = "gbm_files/gbm_stats_healthy.RDS")
saveRDS(object = gmp_fold_rel, file = "gbm_files/gbm_best_tree_healthy.RDS")

# Plot gbm ----------------------------------------------------------------

cm_fold <- readRDS(file = "gbm_files/gbm_stats_healthy.RDS")
gmp_fold_rel <- readRDS(file = "gbm_files/gbm_best_tree_healthy.RDS")

summary(cm_fold)

gmp_fold_rel_means <- apply(X = gmp_fold_rel, MARGIN = 1, function(x) mean(x, na.rm = T))
gmp_fold_rel_medians <- apply(X = gmp_fold_rel, MARGIN = 1, function(x) median(x, na.rm = T))

which(gmp_fold_rel_means  >= 1 ) %>% length
which(gmp_fold_rel_medians  >= 1 ) %>% length

p_gbm_infl <- gmp_fold_rel[names(sort(gmp_fold_rel_means, decreasing = F)), ] %>% 
    rownames_to_column(var = "var") %>% 
    dplyr::filter(var %in% names(which(gmp_fold_rel_means  >= 2 ))) %>% 
    dplyr::mutate(var = gsub(pattern = "B_", replacement = "B: ", x = var)) %>% 
    dplyr::mutate(var = gsub(pattern = "F_", replacement = "F: ", x = var)) %>% 
    dplyr::mutate(var = gsub(pattern = "M_", replacement = "M: ", x = var)) %>% 
    dplyr::mutate(var = gsub(pattern = "__", replacement = " ", x = var)) %>% 
    dplyr::mutate(var = gsub(pattern = "_", replacement = " ", x = var)) %>% 
    dplyr::mutate(var = gsub(pattern = "f x", replacement = "f ", x = var)) %>% 
    dplyr::mutate(var = gsub(pattern = "x uncl", replacement = " uncl", x = var)) %>% 
    pivot_longer(data = ., values_to = "Percentage", cols = starts_with("gbm")) %>% 
    dplyr::filter(!is.na(Percentage)) %>% 
    ggviolin(data = ., x = "var", y = "Percentage", 
             trim = TRUE, fill = "cornflowerblue", alpha = 0.25,
             width = 1.0, orientation = "horiz",
             add = "boxplot", add.params = list(color="grey35", size=0.5)) +
    ylim(0, 100) +
    labs(x = "", y = "Relative influence (%)") +
    ggplot2::annotate(geom = "text", x = 4, y = 75,  
                      label = paste0("Accuracy: ",
                                     round( mean(cm_fold$Accuracy, na.rm = T), 2 ) ), 
                      hjust = 0, parse = TRUE) +
    ggplot2::annotate(geom = "text", x = 3, y = 75,  
                      label = paste0("kappa: ",
                                     round( mean(cm_fold$Kappa, na.rm = T), 2 ) ), 
                      hjust = 0, parse = TRUE) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=12))
p_gbm_infl

mean(cm_fold$Accuracy, na.rm = T); sd(cm_fold$Accuracy, na.rm = T)
mean(cm_fold$Kappa, na.rm = T); sd(cm_fold$Kappa, na.rm = T)

p_gbm_metrics <- cm_fold %>% 
    dplyr::select(Accuracy, Kappa) %>% 
    pivot_longer(data = ., values_to = "value", cols = everything()) %>% 
    dplyr::filter( !is.na(value)) %>% 
    dplyr::filter( value > 0) %>% 
    ggviolin(data = ., x = "name", y = "value", 
             trim = TRUE, fill = "cornflowerblue", alpha = 0.25,
             width = 0.7, orientation = "vertical",
             add = "mean_se") +
    ylim(0, 1.) +
    labs(x = "", y = "") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=12))
    

layout <- "
AAAB
AAAB
AAAB
AAAB
"
p_gbm_infl + p_gbm_metrics +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    plot_annotation(title = 'Healthy: Gradient Boosting Machine (1,000 iterations)') 
