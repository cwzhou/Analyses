library(ggplot2)
library(dplyr)

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
  theme_minimal()

ggplot() +
  # Individual simulation points with a slight jitter for better visualization
  geom_jitter(data = num_re_sim, 
              aes(x = factor(Number_RE), y = count, color = method, group = simulation), 
              alpha = 0.8, size = 1, width = 0.2, height = 0) +
  # Mean points for each Number_RE and method
  geom_point(data = mean_data, 
             aes(x = factor(Number_RE), y = mean_count, color = method), 
             size = 1, shape = 2, stroke = 1.2) +  
  # Dashed line for mean counts across simulations
  geom_line(data = mean_data, 
            aes(x = factor(Number_RE), y = mean_count, color = method, group = method), 
            size = 0.4, linetype = "dashed") +  
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
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("purple", "blue", "red"))

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
                              "zom" = "Zero-Order Model"))  # Custom x-axis labels


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

ggplot(num_re_sim, aes(x = count, 
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
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal"  # Ensure the legend items are laid out horizontally
    ) +
  facet_grid(Number_RE ~ method, scales = "free_y",
             labeller = labeller(
               Number_RE = c("0" = "Total RE: 0",
                             "1" = "Total RE: 1",
                             "2" = "Total RE: 2",
                             "3" = "Total RE: 3"),
               method = c("czmk" = "itrSurv", 
                          "observed" = "observed policy", 
                          "zom" = "zero-order")
             ))  # Change facet labels


# Perform the pivot_longer operations and select only the necessary columns
library(tidyr)
library(dplyr)

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
         survival = 365*survival); result_surv
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


library(ggplot2)
library(patchwork)

# Set factor levels for method
result_long <- result_long %>%
  mutate(method = factor(method, levels = c("itrSurv", "zero-order", "observed")))

# Create the survival jitter plot with means as a line
result_surv_plot <- ggplot(result_long, 
                           aes(x = method, 
                               y = survival,
                               color = factor(sim))) +
  geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jitter for points
  stat_summary(fun = mean, geom = "point", aes(group = 1), 
               color = "black", size = 2) +  # Add means as a line
  labs(x = "Method", y = "Mean Survival Time", 
       title = "Failure Time Comparison", color = "Simulation") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create the MFF jitter plot with means as a line
result_mff_plot <- ggplot(result_long, 
                          aes(x = method, 
                              y = mff,
                              color = factor(sim))) +
  geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jitter for points
  stat_summary(fun = mean, geom = "point", 
               aes(group = 1), color = "black", size = 2) +  # Add means as a line
  labs(x = "Method", y = "MFF", 
       title = "MFF Comparison", color = "Simulation") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine the plots and collect the legend
combined_plot <- (result_surv_plot + result_mff_plot) +
  plot_layout(ncol = 2, guides = "collect") &  # Collect the legends
  theme(legend.position = "bottom")

# Display the combined plot
combined_plot

