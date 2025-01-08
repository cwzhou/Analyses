library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

local = 1
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")
}
date_result = "2025-01-05"
mff_allsims = read.csv(sprintf("output/%s/mff/mff_allsims.csv", date_result))
num_re = read.csv(sprintf("output/%s/mff/num_re.csv", date_result))
result = readRDS(sprintf("output/%s/simResult_RE_censor1_prop1_n1_betaD.1_gammaD.1_omegaD.1_lambda0D.1_tmp.rds", date_result))
mff_allsims = mff_allsims %>% 
  mutate(RE_per_YrsLived = Number_RE/survival)
mff_allsims$method <- factor(mff_allsims$method, levels = c("czmk", "zom", "observed"))

##################################################################################
##################################################################################
##################################################################################
num_re_sim = mff_allsims %>%
  group_by(simulation, Number_RE, method) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(method = factor(method, levels = c("czmk", "zom", "observed")))
head(num_re_sim)
ggplot(num_re_sim, aes(x = factor(Number_RE),
                       y = count,
                       color = method,
                       group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(
    title = "Total Counts of Recurrent Events",
    x = "Total Number of Recurrent Events for a Person",
    y = "Number of People",
    color = "Method"
  ) +
  theme_minimal()

num_re_sim %>% 
  group_by(Number_RE,method) %>% 
  summarise(mean = mean(count), sd = sd(count))

ggplot(num_re_sim, aes(x = as.factor(Number_RE), y = count, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Counts by Number_RE and Method",
    x = "Number_RE",
    y = "Count",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
  ) +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "Observed Policy",
                              "zom" = "Zero-Order Model"))  +
  scale_fill_manual(values = c("purple", "lightblue", "lightgreen"))


# Calculate the mean count for each Number_RE and method
mean_data <- num_re_sim %>%
  group_by(Number_RE, method) %>%
  summarise(mean_count = mean(count), .groups = "drop")

ggplot() +
  # Points for individual simulations
  geom_point(data = num_re_sim,
             aes(x = factor(Number_RE), y = count, color = method, group = simulation),
             alpha = 0.8, size = 2, position = position_jitter(width = 0.2, height = 0)) +
  labs(
    title = "Counts of Total Number of Recurrent Events",
    x = "Total Recurrent Events",
    y = "Frequency",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
  ) +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "Observed Policy",
                              "zom" = "Zero-Order Model"))  +
  scale_color_manual(values = c("purple", "lightblue", "lightgreen"))


ggplot() +
  # Individual simulation points with a slight jitter for better visualization
  geom_jitter(data = num_re_sim,
              aes(x = factor(Number_RE), y = count, color = method, group = simulation),
              alpha = 0.3, size = 1, width = 0.2, height = 0) +
  # Mean points for each Number_RE and method
  geom_point(data = mean_data,
             aes(x = factor(Number_RE), y = mean_count, color = method),
             size = 1, shape = 2, stroke = 1.2) +
  # Dashed line for mean counts across simulations
  geom_line(data = mean_data,
            aes(x = factor(Number_RE), y = mean_count, color = method, group = method),
            size = 1, linetype = "solid") +
  labs(
    title = "Counts of Total Number of Recurrent Events",
    x = "Total Recurrent Events",
    y = "Frequency",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
  ) +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "Observed Policy",
                              "zom" = "Zero-Order Model"))  +
  scale_color_manual(values = c("purple", "lightblue", "lightgreen"))

ggplot(num_re_sim, aes(x = method,
                       y = count,
                       color = factor(Number_RE))) +
  geom_jitter(size = 1,
              width = 0.5,
              shape = 19,
              stroke = 1.2) +  # Add stroke to points for the border color
  labs(
    title = "Count of Total Recurrent Events",
    x = "Method",
    y = "Frequency",
    color = "Total Number of Recurrent Events"  # Renaming the legend title
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
  ) +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "Observed Policy",
                              "zom" = "Zero-Order Model"))  # Custom x-axis labels

ggplot(num_re_sim, aes(x = simulation,
                       y = count,
                       color = factor(simulation))) +
  geom_jitter(size = 0.5,
              width = 0.5,
              shape = 19,
              stroke = 1.2) +  # Add stroke to points for the border color
  labs(
    title = "Total Number of REs Across Simulations",
    x = "",
    y = "Number of People"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.line.x = element_blank(),  # Remove x-axis line
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.major.y = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.y = element_blank(),  # Remove vertical grid lines
    legend.position = "none"
  ) +
  facet_grid(Number_RE ~ method, scales = "free_y",
             labeller = labeller(
               # Number_RE = c("0" = "Total RE: 0",
               #               "1" = "Total RE: 1",
               #               "2" = "Total RE: 2",
               #               "3" = "Total RE: 3"),
               method = c("czmk" = "itrSurv",
                          "observed" = "observed policy",
                          "zom" = "zero-order")
             ))  # Change facet labels

##################################################################################
##################################################################################
##################################################################################
res_sim = mff_allsims %>%
  group_by(simulation, Number_RE, method) %>%
  summarize(count = n(),
            mean_surv = mean(survival),
            mean_RE_yrs = mean(RE_per_YrsLived),
            .groups = "drop") %>%
  mutate(method = factor(method, levels = c("czmk", "zom", "observed")))

ggplot(res_sim, aes(x = as.factor(Number_RE), y = mean_surv, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Distribution of Time to Death",
    x = "Number of Total Recurrent Events for an Individual",
    y = "Survival Time (Years)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
  ) +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "Observed Policy",
                              "zom" = "Zero-Order Model"))  +
  scale_fill_manual(values = c("purple", "lightblue", "lightgreen"))

##################################################################################
##################################################################################
##################################################################################
mean_mff_re = mff_allsims %>%
  group_by(Number_RE, method) %>%
  summarise(mean_survival = mean(survival, na.rm = TRUE),
            mean_RE_yrs = mean(RE_per_YrsLived, na.rm = TRUE))

# Left Plot: Mean Survival vs Method with Number_RE as Facet
plot_survival <- ggplot(mean_mff_re, aes(x = method, y = mean_survival, color = method)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Number_RE, nrow = 1, scales = "free_y") +  # Facet by Number_RE
  scale_y_continuous(
    breaks = seq(0, ceiling(max(mean_mff_re$mean_survival, na.rm = TRUE)), by = 0.1),  # All ticks labeled
    labels = seq(0, ceiling(max(mean_mff_re$mean_survival, na.rm = TRUE)), by = 0.1),
    expand = c(0, 0.05)  # Add padding to ensure grid lines at top and bottom
  ) +
  labs(
    title = "Mean Truncated Survival",
    x = "Method",
    y = "Mean Survival (Years)"
  ) +
  scale_color_manual(
    values = c("czmk" = "#6a0dad", "observed" = "#2ca02c", "zom" = "#1f77b4"),
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "observed policy")
  ) +
  scale_x_discrete(
    limits = c("czmk", "zom", "observed"),  # Order changed here
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "observed policy")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text for readability
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Right Plot: Mean RE per Year vs Method with Number_RE as Facet
plot_re <- ggplot(mean_mff_re, aes(x = method, y = mean_RE_yrs, color = method)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Number_RE, nrow = 1, scales = "free_y") +  # Facet by Number_RE
  scale_y_continuous(
    breaks = seq(0, ceiling(max(mean_mff_re$mean_RE_yrs, na.rm = TRUE)), by = 0.1),  # All ticks labeled
    labels = seq(0, ceiling(max(mean_mff_re$mean_RE_yrs, na.rm = TRUE)), by = 0.1),
    expand = c(0, 0.05)  # Add padding to ensure grid lines at top and bottom
  ) +
  labs(
    title = "Mean Number RE per Years Lived",
    x = "Method",
    y = "Mean Number RE per Years Lived"
  ) +
  scale_color_manual(
    values = c("czmk" = "#6a0dad", "observed" = "#2ca02c", "zom" = "#1f77b4"),
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "observed policy")
  ) +
  scale_x_discrete(
    limits = c("czmk", "zom", "observed"),  # Order changed here
    labels = c("czmk" = "itrSurv", "zom" = "zero-order", "observed" = "observed policy")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text for readability
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Combine the two plots side-by-side with aligned legends
combined_plot <- plot_survival + plot_re +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)

##################################################################################
##################################################################################
##################################################################################
grayscale1 <- c("0" = "white",
                   "1" = "lightgray",
                   "2" = "slategray",  # More distinct contrast with '1'
                   "3" = "black")
# # Plot mean counts for each method, stacked by Number_RE categories
# ggplot(mean_data, aes(x = method, y = mean_count, fill = factor(Number_RE))) +
#   geom_bar(stat = "identity", position = "dodge", color = "black") +  # Black border around bars
#   labs(
#     title = "Mean Count of Total Recurrent Events",
#     x = "Method",
#     y = "Mean Across Simulations",
#     fill = "Total Number of Recurrent Events"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     plot.title = element_text(size = 14, face = "bold"),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 12),
#     panel.grid.major.x = element_blank(),
#     legend.position = "bottom",  # Move legend to the bottom
#     legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
#   ) +
#   scale_fill_manual(values = grayscale1) +
#   scale_x_discrete(labels = c("czmk" = "itrSurv",
#                               "observed" = "observed",
#                               "zom" = "zero-order model"))  # Custom x-axis labels

ggplot(mean_data, aes(x = method,
                      y = mean_count,
                      color = factor(Number_RE))) +
  geom_point(size = 1,
             shape = 19,
             stroke = 1.2) +  # Add stroke to points for the border color
  labs(
    title = "Mean Count of Total Recurrent Events by Method",
    x = "Method",
    y = "Mean Across Simulations",
    color = "Total Number of Recurrent Events"  # Renaming the legend title
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
  ) +
  # scale_fill_manual(values = grayscale1) +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "Observed Policy",
                              "zom" = "Zero-Order Model"))




# Perform the pivot_longer operations and select only the necessary columns
# Reshape and rename levels in the 'method' column
result_surv <- result %>%
  pivot_longer(
    cols = c(czmk_survival, zom_survival, obs_survival),
    names_to = "method",
    values_to = "survival"
  ) %>%
  dplyr::select(sim, method, survival) %>%
  mutate(
    method = case_when(
      method == "czmk_survival" ~ "itrSurv",
      method == "zom_survival" ~ "zero-order",
      method == "obs_survival" ~ "observed",
      TRUE ~ method
    )
  ) %>%
  mutate(method = factor(method, levels = c("itrSurv", "zero-order", "observed")),
         survival = survival); result_surv
result_mff <- result %>%
  pivot_longer(
    cols = c(czmk_endpoint, zom_endpoint, obs_endpoint),
    names_to = "method",
    values_to = "mff"
  ) %>%
  dplyr::select(sim, method, mff) %>%
  mutate(
    method = case_when(
      method == "czmk_endpoint" ~ "itrSurv",
      method == "zom_endpoint" ~ "zero-order",
      method == "obs_endpoint" ~ "observed",
      TRUE ~ method
    )
  ); result_mff
result_long = left_join(result_surv, result_mff, by = c("sim", "method"))
result_long

# Set factor levels for method
result_long <- result_long %>%
  mutate(method = factor(method, levels = c("itrSurv", "zero-order", "observed")))

result_surv_plot <- ggplot(result_long, aes(x = method, y = survival)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jitter for points
  stat_summary(fun = mean,
               geom = "point",
               aes(color = method, group = 1),  # Map color to method
               shape = "square",
               size = 3) +  # Add means as points
  stat_summary(fun = mean, geom = "text",
               aes(color = method, label = round(after_stat(y), 2)),
               vjust = 1.5, size = 4) +  # Display mean values as text
  labs(x = "Method", y = "Restricted Mean Survival Time at Tau",
       title = "Phase 1 Endpoint Comparison", color = "Simulation") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("czmk" = "itrSurv", "observed" = "Observed Policy", "zom" = "Zero-Order Model"))
result_surv_plot


# Create the MFF jitter plot with means as a line
result_mff_plot <- ggplot(result_long,
                          aes(x = method,
                              y = mff)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jitter for points
  stat_summary(fun = mean,
               geom = "point",
               aes(color = method, group = 1),  # Map color to method
               shape = "square",
               size = 3) +  # Add means as points
  stat_summary(fun = mean, geom = "text",
               aes(color = method,
                   label = round(after_stat(y), 2)),
               vjust = 1.5, size = 4) +
  labs(x = "Method", y = "Average Mean Frequency at Tau",
       title = "Phase 2 Endpoint Comparison") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("czmk" = "itrSurv",
                              "observed" = "observed policy",
                              "zom" = "zero-order"))
result_mff_plot

# Combine the plots and collect the legend
combined_plot <- (result_surv_plot + result_mff_plot) +
  plot_layout(ncol = 2, guides = "collect") &  # Collect the legends
  theme(legend.position = "none")
# Display the combined plot
combined_plot


##################################################################################
##################################################################################
##################################################################################
##################################################################################
library(dplyr)
library(purrr)
row_re_datasets <- function(prefix = "rep", suffix = "sim1") {
  # Dynamically construct dataset names
  datasets <- list(
    rep_czmk = get(paste0(prefix, "_czmk_", suffix))$dataset_recurrent,
    rep_zom = get(paste0(prefix, "_zom_", suffix))$dataset_recurrent,
    rep_obs = get(paste0(prefix, "_obs_", suffix))$dataset_recurrent
    )
  # Calculate row counts for all recurrent levels in each dataset
  row_counts <- map_dfr(datasets, function(dataset) {
    dataset %>%
      group_by(R_Label) %>%
      summarise(RowCount = n(), .groups = "drop")
  }, .id = "Dataset")

  return(row_counts)
}

# Initialize lists to store results
results_long <- list()
results_wide <- list()
# Process each simulation
for (i in 1:n.sim) {
  suffix <- paste0("sim", i)

  # Process datasets for the current simulation
  result_sim <- row_re_datasets(prefix = "rep", suffix = suffix)
  results_long[[i]] <- result_sim

  # Convert to wide format
  result_sim_wide <- result_sim %>%
    pivot_wider(
      names_from = Dataset,  # Create columns for each dataset
      values_from = RowCount # Fill values from the RowCount column
    )
  results_wide[[i]] <- result_sim_wide
}
print(results_long[[1]])  # Long format for sim1
print(results_wide[[1]])  # Wide format for sim1
print(results_long[[2]])  # Long format for sim1
print(results_wide[[2]])  # Wide format for sim1
