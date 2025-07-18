library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

read_data = 0
local = 1
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")
}
# date_result_list = c("2025-01-14", "2025-01-15")
# for (date_result in date_result_list){
#   mff_allsims = read.csv(sprintf("output/%s/mff/mff_allsims.csv", date_result))
#   # num_re = read.csv(sprintf("output/%s/mff/num_re.csv", date_result))
#   # result = readRDS(sprintf("output/%s/simResult_RE_censor1_prop2_n1_betaD.1_gammaD.1_omegaD.1_lambda0D.1_tmp.rds", date_result))
# }

if (read_rds == 1){
  home_dir = "~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE"
  date_result_list <- "2025-03-05"
  file_path <- file.path(home_dir, "output", date_result_list, "mff")
  
  # Generate file indices, excluding 801-900
  file_indices <- setdiff(1:1100, 801:900)
  
  # Generate full file paths
  file_names <- file.path(file_path, paste0("stacked.mff_sim", file_indices, ".rds"))
  
  # Read all .rds files into a list
  data_list <- lapply(file_names, readRDS)
  
  # Optionally combine if they have the same structure
  library(dplyr)
  mff_allsims0 <- bind_rows(data_list)
  # View(mff_allsims0)
  mff_allsims = mff_allsims0 %>%
    # mutate(RE_cat = ifelse(Number_RE <= 1.5, 'Low',
    #                        ifelse(Number_RE > 2.5 & Number_RE <= 4, 'Med',
    #                               'High'))) %>%
    # mutate(RE_cat = factor(RE_cat, levels = c("Low", "Med", "High"))) %>%
    mutate(RE_cat = ifelse(Number_RE <= 2, '[0,2]','(2,4]')) %>%
    mutate(RE_cat = factor(RE_cat, levels = c('[0,2]','(2,4]'))) %>%
    mutate(RE_per_YrsLived = Number_RE/survival)
  mff_allsims$method <- factor(mff_allsims$method, levels = c("czmk", "zom", "observed"))
  
}
if (read_data == 1){
  date_result_list <- "2025-03-05"
    # c("2025-03-01",
    #   "2025-03-02","2025-03-03",
    #   "2025-03-04","2025-03-05","2025-03-06",
    #   "2025-03-07","2025-03-08","2025-03-09",
    #   "2025-03-10","2025-03-11")
    # c("2025-03-02.1", "2025-03-02.2", "2025-03-02.3", "2025-03-02.4", "2025-03-02.5", 
    #   "2025-03-02.6", "2025-03-02.7", "2025-03-02.8", "2025-03-02.9"
    #                     )
    # c("2025-02-01", "2025-02-02", "2025-02-03", "2025-02-04",
    #                     "2025-02-05", "2025-02-06", "2025-02-07", "2025-02-09",
    #                     "2025-02-10", "2025-02-11")
    # c("2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04", "2025-01-05", "2025-01-06",
    # "2025-01-07", #"2025-01-18",
    # "2025-01-09", "2025-01-10", "2025-01-11")
  # c(#"2025-01-11", 
  #   "2025-01-13", "2025-01-14", "2025-01-15", "2025-01-16",
  # "2025-01-17", "2025-01-18", #"2025-01-19",
  # "2025-01-20", "2025-01-21") # need to add "older" after output/
  mff_sims_list <- list()
  for (date_result in date_result_list) {
    file_path <- sprintf("output/%s/mff/mff_allsims.csv", date_result)
    if (file.exists(file_path)) {
      mff_sims <- read.csv(file_path)
      mff_sims_list[[length(mff_sims_list) + 1]] <- mff_sims
    } else {
      warning(sprintf("File not found: %s", file_path))
    }
  }
  # Combine all data frames into a single data frame
  mff_allsims0 <- do.call(rbind, mff_sims_list)
  # print(head(mff_allsims0)); print(tail(mff_allsims0))
  # print(length(unique(mff_allsims0$simulation)))
  mff_allsims = mff_allsims0 %>%
    # mutate(RE_cat = ifelse(Number_RE <= 1.5, 'Low',
    #                        ifelse(Number_RE > 2.5 & Number_RE <= 4, 'Med',
    #                               'High'))) %>%
    # mutate(RE_cat = factor(RE_cat, levels = c("Low", "Med", "High"))) %>%
    mutate(RE_cat = ifelse(Number_RE <= 2, '[0,2]','(2,4]')) %>%
    mutate(RE_cat = factor(RE_cat, levels = c('[0,2]','(2,4]'))) %>%
    mutate(RE_per_YrsLived = Number_RE/survival)
  mff_allsims$method <- factor(mff_allsims$method, levels = c("czmk", "zom", "observed"))
} else{
  message("not reading data because already done")
}


n.sim = length(unique(mff_allsims0$simulation))
mff_00 = mff_allsims

mean1 = mff_allsims %>%
  group_by(simulation, method) %>%
  summarize(meanRE = mean(Number_RE),
            meanSURV = mean(survival),
            meanTRT1 = mean(Trt),
            meanREYRS = mean(RE_lived))

mean2 = mean1 %>%
  group_by(method) %>%
  summarise(meanRE1 = mean(meanRE),
            meanSURV1 = mean(meanSURV),
            meanTRT11 = mean(meanTRT1),
            meanREYRS1 = mean(meanREYRS))

# View(mean1)
# View(mean2)

mff_allsims0 %>%
  group_by(method) %>%
  summarise(
    #n = n(),
    meanSurv = mean(survival),
    sdSurv = sd(survival),
    Q1 = quantile(survival, 0.25),
    Q2 = quantile(survival, 0.5),  # Median (Q2)
    Q3 = quantile(survival, 0.75)
  )

mff_allsims0 %>%
  group_by(method) %>%
  summarise(
    #n = n(),
    meanREyr = mean(RE_lived),
    sdREyr = sd(RE_lived),
    Q1 = quantile(RE_lived, 0.25),
    Q2 = quantile(RE_lived, 0.5),  # Median (Q2)
    Q3 = quantile(RE_lived, 0.75)
  )


mff_allsims0 %>%
  group_by(method) %>%
  summarise(
    #n = n(),
    meanTrt1 = mean(Trt),
    meanSurv = mean(survival),
    sdTrt = sd(Trt),
    Q1 = quantile(Trt, 0.25),
    Q2 = quantile(Trt, 0.5),  # Median (Q2)
    Q3 = quantile(Trt, 0.75)
  )

mff_allsims0 %>%
  group_by(method, Number_RE) %>%
  summarise(
    #n = n(),
    meanTrt = mean(Trt)
  )

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
# Summary statistics
mff_means = mff_allsims %>%
  group_by(simulation, method) %>%
  summarise(mean_surv = mean(survival),
            mean_re_yrs = mean(RE_per_YrsLived))
# View(mff_means)

# Get global y-axis limits for both plots
ymax = max(mff_means$mean_re_yrs, na.rm = TRUE)
y_max <- max(ceiling(max(mff_means$mean_surv, na.rm = TRUE)),
             ymax)
y_min <- 0.5  # Start from 0 for both plots
set.seed(1)
# Left Plot: Mean Survival vs Method with Number_RE as Facet
plot_survival <- ggplot(mff_means, aes(x = method, y = mean_surv, color = method)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, height = 0, size = 1) +  # Add jitter for points
  stat_summary(fun = mean,
               geom = "point",
               aes(group = 1),  # Map color to method
               shape = "square", size = 3, color = "black") +
  stat_summary(fun = mean, geom = "text",
               aes(label = round(after_stat(y), 2)),
               vjust = -2, size = 5, color = "black") +
  # scale_y_continuous(
  #   limits = c(y_min, y_max),  # Set the same y-axis limits for both plots
  #   breaks = seq(y_min, y_max, by = 0.5),
  #   labels = seq(y_min, y_max, by = 0.5),
  #   expand = c(0, 0.05)
  # ) +
  labs(
    title = "Phase 1 Endpoint Comparison",
    x = "",
    y = "Truncated Mean Survival (Years)"
  ) +
  scale_color_manual(
    # values = c("czmk" = "#6a0dad", "zom" = "#1f77b4", "observed" = "#2ca02c"),
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  ) +
  scale_x_discrete(
    limits = c("czmk", "zom", "observed"),  # Order changed here
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 15, hjust = 1)  # Rotate x-axis text for readability
  )

# Right Plot: Mean RE per Year vs Method with Number_RE as Facet
plot_re <- ggplot(mff_means, aes(x = method, y = mean_re_yrs, color = method)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, height = 0, size = 1) +  # Add jitter for points
  stat_summary(fun = mean,
               geom = "point",
               aes(group = 1),  # Map color to method
               shape = "square", size = 3, color = "black") +
  stat_summary(fun = mean, geom = "text",
               aes(label = round(after_stat(y), 2)),
               vjust = -2, size = 5, color = "black") +
  # scale_y_continuous(
  #   limits = c(y_min, y_max),  # Set the same y-axis limits for both plots
  #   breaks = seq(y_min, y_max, by = 0.5),
  #   labels = seq(y_min, y_max, by = 0.5),
  #   expand = c(0, 0.05)
  # ) +
  labs(
    title = "Phase 2 Endpoint Comparison",
    x = "",
    y = "Mean Number RE\nper Year at Risk"
  ) +
  scale_color_manual(
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  ) +
  scale_x_discrete(
    limits = c("czmk", "zom", "observed"),  # Order changed here
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 15, hjust = 1)  # Rotate x-axis text for readability
  )

# Combine the two plots side-by-side with aligned legends
combined_plot <- plot_survival + plot_re +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom"); print(combined_plot)


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
# by Number_RE
mff_means1 = mff_allsims %>%
  group_by(simulation, Number_RE, method) %>%
  summarise(mean_surv = mean(survival),
            mean_re_yrs = mean(RE_per_YrsLived))
# mff_means1 = mff_allsims %>%
#   group_by(simulation, RE_cat, method) %>%
#     summarise(mean_surv = mean(survival),
#               mean_re_yrs = mean(RE_per_YrsLived))
# mff_means1 = mff_allsims %>%
#   group_by(method, RE_cat) %>%
#   summarise(mean_surv = mean(survival),
#             mean_re_yrs = mean(RE_per_YrsLived),
#             n = n())

# head(mff_means1)
set.seed(11)
# Left Plot: Mean Survival vs Method with Number_RE as Facet
plot_survival1 <-  ggplot(mff_means1,
                          aes(x = factor(Number_RE), #as.factor(RE_cat),
                              y = mean_surv,
                              color = method)) +
  geom_boxplot(position = position_dodge(0.8)) +  # Adjust the dodge width if necessary
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), size = 1) +  # Adjust jitter dodge
  stat_summary(fun = mean, geom = "point",
               aes(group = method),
               position = position_dodge(0.8),
               color = "black", size = 3, shape = 18) +  # Adds mean points per method
  stat_summary(fun = mean, geom = "text",
               aes(label = sprintf("%.3f", ..y..),
                   group = method),
               position = position_dodge(0.8),
               color = "black", size = 3, vjust = -1) +
  labs(
    title = "Phase 1 Endpoint Comparison",
    x = "Total Number of Recurrent Events",
    y = "Mean Truncated\nSurvival in years"
  ) +
  scale_color_manual(
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  ) +
  # scale_x_discrete(
  #   limits = c("czmk", "zom", "observed"),  # Order changed here
  #   labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  # ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_text(angle = 15, hjust = 1)  # Rotate x-axis text for readability
  )

# Right Plot: Mean RE per Year vs Method with Number_RE as Facet
plot_re1 <- ggplot(mff_means1,
                   aes(x = factor(Number_RE),
                       y = mean_re_yrs,
                       color = method)) +
  geom_boxplot(position = position_dodge(0.8)) +  # Adjust the dodge width if necessary
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), size = 1) +  # Adjust jitter dodge
  stat_summary(fun = mean, geom = "point",
               aes(group = method),
               position = position_dodge(0.8),
               color = "black", size = 3, shape = 18) +  # Adds mean points per method
  stat_summary(fun = mean, geom = "text",
               aes(label = sprintf("%.3f", ..y..),
                   group = method),
               position = position_dodge(0.8),
               color = "black", size = 3, vjust = -1) +
  labs(
    title = "Phase 2 Endpoint Comparison",
    x = "Total Number of Recurrent Events",
    y = "Mean Number RE\nper Year at Risk"
  ) +
  scale_color_manual(
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  ) +
  # scale_x_discrete(
  #   limits = c("czmk", "zom", "observed"),  # Order changed here
  #   labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  # ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_text(angle = 15, hjust = 1)  # Rotate x-axis text for readability
  )

# Combine the two plots side-by-side with aligned legends
combined_plot1 <- plot_survival1 + plot_re1 +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom"); print(combined_plot1)



#############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
summary = mff_allsims %>%
  mutate(RE_cat = ifelse(Number_RE <= 2, 'Low',
                         ifelse(Number_RE > 2 & Number_RE <= 6, 'Med',
                                'High'))) %>%
  mutate(RE_cat = factor(RE_cat, levels = c("Low", "Med", "High"))) %>%
  group_by(method, RE_cat) %>% #Number_RE) %>%
  summarise(#min = min(RE_per_YrsLived),
            #max = max(RE_per_YrsLived),
            mean_RE_yrs = mean(RE_per_YrsLived),
            sd_RE_yrs = sd(RE_per_YrsLived),
            mean_surv = mean(survival))
# View(summary)

ps = ggplot(summary, aes(
  x = factor(RE_cat),
  y = mean_surv,
  color = method,
  group = method
)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +  # Dodge points
  # scale_color_brewer(palette = "Set2") +  # Improved color palette
  labs(
    title = "",
    x = "Total Number of RE",
    y = "Truncated Mean Survival (years)"
  ) +
  theme_minimal(base_size = 14) +  # Larger base font size for readability
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  )

pre = ggplot(summary, aes(
  x = factor(RE_cat),
  y = mean_RE_yrs,
  color = method,
  group = method
)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +  # Dodge points
  # scale_color_brewer(palette = "Set2") +  # Improved color palette
  labs(
    title = "",
    x = "Total Number of RE",
    y = "Mean RE\nper Year at Risk"
  ) +
  theme_minimal(base_size = 14) +  # Larger base font size for readability
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  )

# Combine the two plots side-by-side with aligned legends
cp <- ps + pre +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom"); print(cp)

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

ggplot(summary, aes(x = method,
                    y = mean_RE_yrs,
                    color = method,
                    group = method)) +
  geom_point(size = 3) +  # Larger points
  # geom_line(size = 1.2) +  # Thicker lines
  # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, size = 1) +  # Wider and thicker error bars
  facet_wrap(~factor(RE_cat), scales = "free_y") +  # Facets by method
  # scale_color_brewer(palette = "Set2") +  # Improved color palette
  labs(
    title = "",
    x = "",
    y = "Mean RE\nper Year at Risk"
  ) +
  theme_minimal(base_size = 14) +  # Larger base font size for readability
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("czmk" = "#c39bd8",  # Lighter purple version of #6a0dad
               "zom" = "#80b0e0",  # Lighter blue
               "observed" = "#7bcf7b"),  # Lighter green
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "randomized policy")
  )


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

