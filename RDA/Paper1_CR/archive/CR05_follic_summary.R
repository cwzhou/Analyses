library(dplyr); library(tibble); library(tidyverse)
library(pracma) # for linking NA values in geom_lines
library(ggplot2); library(cowplot)

# add in other comparator methods

  # Specify date of the outputs and all_cause (T/F)
  date = "2024-01-27" #Sys.Date()
  
  mthd = c("CZMK", "CSK", "PMCR", "ZOM", "CSKzom", "observed") 
  labs = c("the proposed method", "Cho et al (2022)", "PMCR (2021)", 
           #"Goldberg & Kosorok (2012), RF", "Goldberg & Kosorok (2012), linear", "Simoneau et al. (2019)", 
           "zero-order model", "Cho zero-order model", "observed")
  cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  names(cols) = labs
  crit = "mean"
  nodesize = 50
  Tx = "A"
  Tx.full = "Treatment"
  clipping = "5%"  # clipping at 5%
  endpoint = "CR"
  dataset_name = "follic"
  K = K
  tau = tau

  p1 = p2 = list()
  for (crit in c("mean", "mean.surv.combo")) {
    # OSValues_20CV_mean_rulemean_tau32_2024-01-26.rds
    nm.tmp =
      sprintf("../3_output/%s/%s/%s/OSvalues_%sCV_%s_tau%s_%s.rds", 
              endpoint, dataset_name, date, K,
              if (crit == "mean") "mean_rulemean" else "surv.mean2200_rulelogrank", 
              tau, date) 
    data.long = 
      nm.tmp %>%
      readRDS() %>% 
      as.data.frame %>% 
      dplyr::select(-train_cens, -train_cause1, -train_cause2, -ns.CZMK, -ns.CSK,) %>% 
      mutate(rep = 1:n()) %>% 
      pivot_longer(cols = -rep,
        names_to = c("method", "var"),
        names_sep = "\\."
      ) %>%
      pivot_wider(
        names_from = var,
        values_from = value
      ) %>%
      # tidyr::gather(key = "method", value = "OS", -rep, -contains(".PC")) %>% 
      # tidyr::gather(key = "method", value = "Cause1", -rep, -method, -OS) %>% 
      mutate(method = factor(method, 
                             levels = mthd,
                             labels = labs),
             crit = crit)
    
    plot_RDA <- function(data, endpoint_type) {
      ggplot(data) +
        geom_boxplot(aes(method, !!sym(endpoint_type), col = method)) +
        geom_jitter(aes(method, !!sym(endpoint_type), col = method, group = rep), 
                    width = 0.2, alpha = 1) +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        scale_color_manual(values = cols) + 
        stat_summary(aes(x = as.numeric(method), 
                         y = !!sym(endpoint_type)),
                     fun.y = mean, 
                     geom = 'point', 
                     col = "black", 
                     shape = "square", 
                     size = 1) +
        stat_summary(aes(x = as.numeric(method) + 0.2, 
                         y = !!sym(endpoint_type), 
                         label = sprintf( paste0("%3.", if (crit == "mean") 0 else 3, "f"), ..y..)),
                     fun.y = mean, 
                     geom = 'text',  
                     col = "black", 
                     position = 'dodge', 
                     vjust = 1) +
        theme_bw() + theme(axis.title.x = element_blank()) +
        # ylim(c(if (crit == "mean") { if (all_cause) 2000 else 2200} else { if (all_cause) 0.45 else 0.6}, NA)) +
        ylab(if (crit == "mean"){
          if (endpoint_type == "OS"){
            "Truncated mean overall survival time"
          } else{
            sprintf("Truncated mean cause %s survival time", priority_cause)
          }
          } else{
            if (endpoint_type == "PC"){
              "Overall truncated survival probability"
            } else{
              sprintf("Overall truncated cause %s survival probability", priority_cause)
            }
            }) +
        guides(col = "none", group = "none")
    }
    
    p1[[crit]] = plot_RDA(data.long, "OS")
    p2[[crit]] = plot_RDA(data.long, "PC")
    
    p.grid = plot_grid(p1[[crit]] + guides(col = FALSE),
              p2[[crit]] + guides(col = FALSE),
              align = "h", nrow = 1, ncol = 2, common.legend = TRUE)
    p.grid2[[crit]] <- plot_grid(p.grid, get_legend(p1[[1]]), align = "v", nrow = 2, rel_heights = c(1, 0.1))
    
    ggsave(nm.tmp %>% 
             gsub("output/", "figure/", .) %>% 
             gsub("OSvalues", "Figure", .) %>% 
             gsub("\\.rds", "\\.png", .),
           width = 10, height = 4)
  }
  