# dev.off();
library(dplyr); library(tibble); library(tidyverse)
library(pracma) # for linking NA values in geom_lines
library(ggplot2); library(cowplot)

# add in other comparator methods

# Specify date of the outputs and all_cause (T/F)
date = "2024-08-25";#"2024-03-03";#"2024-02-26"; #"2024-01-28" #Sys.Date()

mthd0 = c("CZMK", "CSK", "PMCR", "AIPWE", "ZOM", "CSKzom", "observed")
labs0 = c("itrSurv", "dtrSurv (2023)", "PMCR (2021)", "AIPWE (2021)",
         "zero-order model", "Cho zero-order model", "observed")
cols0 <- c("#614CFF", "#00BFC9", "#F8766D", "#00BA38", "#FFA500", "#F564E3", "#800040")
names(cols0) = labs0

keep_method <- c(TRUE,TRUE,TRUE,TRUE,TRUE,!TRUE,TRUE);
mthd <- mthd0[keep_method]
labs <- labs0[keep_method]
cols <- cols0[keep_method]
# Concatenate ".OS" to each element of the vector
mthd.OS <- paste(mthd, "OS", sep = ".")
# Concatenate ".PC" to each element of the vector
mthd.PC <- paste(mthd, "PC", sep = ".")
# Combine the two vectors into endpoint vector
mthd.ep <- c(mthd.OS, mthd.PC)

possible_crits1 = c("area", "mean", "mean.prob.combo", "prob")
possible_crits1 = possible_crits1[2]
# crit1 = 1

nodesize = 5#50
Tx = "A"
Tx.full = "Treatment"
clipping = "5%"  # clipping at 5%
endpoint = "CR"
dataset_name = "pad"
K = K
tau = tau
p1 = p2 = list()


for (crit1 in 1:length(possible_crits1)) {
  crit = possible_crits[crit1]
  # Values_200CV_mean_rule1mean_surv_rule2gray_cr_tau365_2024-08-25.rds
  # for now, this is fixed rule2 (gray's test)
  # for now, this is fixed CR
  nm.tmp =
    sprintf("./3_output/%s/%s/%s/Values_%sCV_%s_rule2gray_cr_tau%s_%s.rds",
            endpoint, dataset_name, date, K,
            if (crit == "mean") "mean_rule1mean_surv"
            else if (crit == "area") "area_rule1area_surv"
            else if (crit == "prob") sprintf("prob%s_rule1logrank_surv",t0_crit)
            else sprintf("mean.prob.combo%s_rule1logrank_surv",t0_crit), #"surv.mean2200_rulelogrank",
            tau, date)
  data.long =
    nm.tmp %>%
    readRDS() %>%
    as.data.frame %>%
    #do not select columns where all rows are NA (aka if we dont run that method)
    dplyr::select(where(~ !all(is.na(.)))) %>%
    # dplyr::select(-train_cens, -train_cause1, -train_cause2, -ns.CZMK, -ns.CSK,) %>%
    dplyr::select(all_of(mthd.ep)) %>%
    mutate(rep = 1:n()) %>%
    pivot_longer(cols = -rep,
      names_to = c("method", "var"),
      names_sep = "\\."
    ) %>%
    pivot_wider(
      names_from = var,
      values_from = value
    ) %>%
    mutate(method = factor(method,
                           levels = mthd,
                           labels = labs),
           crit = crit)
  # note: AIPWE and PMCR have some NaN; AIPWE is missing sometimes

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
                   size = 2) +
      # stat_summary(aes(x = as.numeric(method),
      #                  y = !!sym(endpoint_type),
      #                  label = sprintf( paste0("%3.", if (crit == "mean") 0 else 3, "f"), ..y..)),
      #              fun.y = mean,
      #              geom = 'text',
      #              col = "black",
      #              position = 'dodge',
      #              vjust = -10,
      #              fontface = "bold",
      #              size = 6) +
      ggplot2::stat_summary(
        aes(
          x = as.numeric(method),
          y = !!sym(endpoint_type),
          label = sprintf(paste0("%3.", if (crit == "mean") 0 else 3, "f"), ..y..)
        ),
        fun.y = mean,
        geom = 'label',  # Use geom_label for a background
        col = "black",
        fill = "white",  # White background
        label.size = 0,  # No border
        position = position_dodge(width = 0.9),  # Adjust the position to avoid overlap
        vjust = -0.5,
        size = 3,
        fontface = "bold",  # Make the text bold
        label.padding = unit(0.25, "lines")  # Padding around the label
      ) +
      theme_bw() + theme(axis.title.x = element_blank()) +
      # ylim(c(if (crit == "mean") { if (all_cause) 2000 else 2200} else { if (all_cause) 0.45 else 0.6}, NA)) +
      ylab(if (crit == "mean" | crit == "area"){
        if (endpoint_type == "OS"){
          "Mean truncated \noverall survival (days)"
        } else{
          "Mean truncated \ndeath survival (days)"
        }
        } else{
          if (endpoint_type == "PC"){
            "Truncated overall survival probability"
          } else{
            "Truncated death survival probability?"
          }
          }) +
      guides(col = "none", group = "none")
  }

  p1[[crit1]] = plot_RDA(data.long, "OS")
  p2[[crit1]] = plot_RDA(data.long, "PC")

  p.grid = plot_grid(p1[[crit1]] + guides(col = FALSE),
            p2[[crit1]] + guides(col = FALSE),
            align = "v",
            nrow = 1, ncol = 2,
            common.legend = TRUE,
            axis = "v")
  # p.grid2 <- plot_grid(p.grid,
  #                      get_legend(p1[[crit1]]),
  #                      align = "v", nrow = 2,
  #                      rel_heights = c(1, 0.1))
  # dont need above b/c no axis to share
  ggsave(nm.tmp %>%
           gsub("output/", "figure/", .) %>%
           gsub("Values", "Figure", .) %>%
           gsub("\\.rds", "\\.png", .),
         width = 10, height = 4)
}
