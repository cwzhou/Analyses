# dev.off();
library(dplyr); library(tibble); library(tidyverse)
library(pracma) # for linking NA values in geom_lines
library(ggplot2); library(cowplot)

mthd0 = c("CZMK", "ZOM", "OBS")
labs0 = c("itrSurv", "zero-order model", "observed policy")
cols0 <- c("#614CFF", "#00BFC9",  "#00BA38")
names(cols0) = labs0

keep_method <- c(TRUE,TRUE,TRUE);
mthd <- mthd0[keep_method]
labs <- labs0[keep_method]
cols <- cols0[keep_method]
# Concatenate ".OS" to each element of the vector
mthd.OS <- paste(mthd, "terminal", sep = ".")
# Concatenate ".PC" to each element of the vector
mthd.PC <- paste(mthd, "RE", sep = ".")
# Combine the two vectors into endpoint vector
mthd.ep <- c(mthd.OS, mthd.PC)

possible_crits1 = c("area", "mean", "mean.prob.combo", "prob")
crit = possible_crits1[2]
crit1 = 1

Tx.full = "Treatment"
clipping = "5%"  # clipping at 5%
endpoint = "RE"
dataset_name = "bladder"
p1 = p2 = list()

data.long = values %>%
  as.data.frame %>%
  #do not select columns where all rows are NA (aka if we dont run that method)
  dplyr::select(where(~ !all(is.na(.)))) %>%
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
      ggplot2::stat_summary(
        aes(
          x = as.numeric(method),
          y = !!sym(endpoint_type),
          label = sprintf(paste0("%3.", if (endpoint_type == "RE") 4 else 1, "f"), ..y..)
        ),
        fun.y = mean,
        geom = 'text',
        col = "black",
        position = position_dodge(width = 0.9),  # Adjust the position to avoid overlap
        vjust = 1.7,
        # hjust = -0.2,
        size = 3,
        fontface = "bold",  # Make the text bold
      ) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            panel.grid.major = element_blank()#, # Removes major gridlines
            # panel.grid.minor = element_blank()  # Removes minor gridlines
      ) +
      # ylim(c(if (crit == "mean") { if (all_cause) 2000 else 2200} else { if (all_cause) 0.45 else 0.6}, NA)) +
      ylab(if (crit == "mean" | crit == "area"){
        if (endpoint_type == "terminal"){
          "Mean truncated death time (months)"
        } else{
          "Mean cancer recurrence\nper truncated months lived"
        }
      } else{
        stop("Not coded up yet for other crit values")
      }) +
      guides(col = "none", group = "none") #+
  }

  p1[[crit1]] = plot_RDA(data.long, "terminal"); p1[[crit1]]
  p2[[crit1]] = plot_RDA(data.long, "RE"); p2[[crit1]]

  p.grid = plot_grid(p1[[crit1]] + guides(col = FALSE),
                     p2[[crit1]] + guides(col = FALSE),
                     align = "v",
                     nrow = 1,
                     ncol = 2,
                     # common.legend = TRUE,
                     axis = "v");

  # Add a title using ggdraw() and draw_label()
  p.final <- ggdraw() +
    draw_label("Figure X. Application to Bladder Dataset using 5-Fold Cross-Validation", x = 0.5, y = 0.95, size = 10, hjust = 0.5) +
    draw_plot(p.grid, x = 0, y = 0, width = 1, height = 0.9);print(p.final)

  ggsave(paste0(fnm %>%
           gsub("output/", "figure/", .), "BladderFigure_", K, "CV.png"),
         width = 10, height = 4)

  