
# libraries ---------------------------------------------------------------

library(tidyverse)
library(KODAMA)

# Data and colors ---------------------------------------------------------

colv <- c("steelblue", "orange")

cova   <- readRDS(file = "data/processed_cova.RDS")

cova$group <- "DM II"
cova$group[grep(pattern = "control", x = cova$Status) ] <- "healthy"

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

# Run KODAMA (Knowledge Discovery by Accuracy Maximization) ---------------

k_data <- cbind(norm_b.clr, norm_f.clr, norm_m.clr)
sum( colSums(k_data) == 0 )
kodama_pls_test <- list()
kodama_pls_test[["2"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 2)
kodama_pls_test[["4"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 4)
kodama_pls_test[["6"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 6)
kodama_pls_test[["8"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 8)
kodama_pls_test[["10"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 10)
kodama_pls_test[["20"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 20)
kodama_pls_test[["50"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 50)
kodama_pls_test[["100"]] <- KODAMA.matrix(data = k_data, FUN = "PLS-DA", f.par = 100)

kodama.table <-  data.frame( 
    Parameter = c(2,4,6,8,10,20,50,100), 
    Entropy = c(kodama_pls_test[["2"]]$entropy, kodama_pls_test[["4"]]$entropy, 
    kodama_pls_test[["6"]]$entropy, kodama_pls_test[["8"]]$entropy,
    kodama_pls_test[["10"]]$entropy, kodama_pls_test[["20"]]$entropy,
    kodama_pls_test[["50"]]$entropy, kodama_pls_test[["100"]]$entropy) ) 
# saveRDS(object = kodama.table, file = "data/KODAMA.entropy.RDS")

kodama.table <- readRDS(file = "data/KODAMA.entropy.RDS")
knitr::kable(kodama.table)

# Minimum Entropy?
kodama.table %>% dplyr::filter(Entropy == min(Entropy))

kodama_pls <- KODAMA.matrix(data = k_data, 
                            FUN = "PLS-DA", 
                            f.par = 8, M = 100, 
                            bagging = TRUE, dims = 4,
                            constrain = factor(cova$group) )

kodama_pls$entropy

kkplot <- KODAMA.visualization(kk = kodama_pls, method = "MDS")


# k.test(data = kkplot, labels = as.numeric(factor(cova$Fiber.Content.Diet)), 
       # n = 10000)
# k.test(data = kkplot, labels = as.numeric(factor(cova$group)), 
       # n = 10000)

k_data_plot <- kkplot %>% 
    data.frame() %>% 
    dplyr::mutate(Diet = cova$Fiber.Content.Diet) %>% 
    dplyr::mutate(Group = cova$group) %>% 
    dplyr::rename(X1 = First.Dimension) %>% 
    dplyr::rename(X2 = Second.Dimension) %>% 
    tibble()

k_data_plot %>% 
    dplyr::mutate(Group = factor(x = Group, levels = c("healthy", "DM II"))) %>% 
    dplyr::mutate(Diet  = gsub(pattern = "High fiber diet", replacement = "HFD", x = Diet)) %>% 
    dplyr::mutate(Diet  = gsub(pattern = "Low fiber diet", replacement = "LFD", x = Diet)) %>% 
    ggplot(data = ., mapping = aes(x = X1, y = X2, color = Diet, 
                                   shape = Group ) ) +
    geom_point() +
    scale_color_manual(values=colv) +
    # scale_color_continuous(type = "viridis") +
    geom_point(size=4) +
    # ggtitle( "Semi-supervised KODAMA with constraint" ) +
    ggtitle( "" ) +
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
ggsave(filename = "plots_final/Figure_4.pdf", height = 6, width = 6)

m_h1 <- lm(formula = X1 ~ Diet, data = k_data_plot %>% dplyr::filter(Group == "healthy"))
m_h2 <- lm(formula = X2 ~ Diet, data = k_data_plot %>% dplyr::filter(Group == "healthy"))
m_d1 <- lm(formula = X1 ~ Diet, data = k_data_plot %>% dplyr::filter(Group == "DM II"))
m_d2 <- lm(formula = X2 ~ Diet, data = k_data_plot %>% dplyr::filter(Group == "DM II"))

summary(m_h1)
summary(m_h2)

summary(m_d1)
summary(m_d2)


data.frame( loadings = loads(kodama_pls), feature = colnames(k_data) ) %>% 
    arrange( desc(loadings) ) %>% 
    head(n = 10)
