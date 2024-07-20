# Define the methods and their corresponding names
methods = method.nm.simple
eval_names <- c("OS_eval", "CIF_eval")
action_names <- "action"
rename_mappings <- c("OS_eval" = "_OS",
                     "CIF_eval" = "_CIF",
                     "action" = "_action")

combined_data_long = list()
for (sim in 1:n.sim){
  sim_lab = sprintf("sim%s", sim)
  # Initialize an empty list to store dataframes
  rep_dfs1 <- list()

  # Loop over each method
  for (method in methods) {
    # Get the corresponding dataframe
    eval_dat = full[["true_P2_eval"]][[sim_lab]][["eval_datasets"]]
    rep_df <- eval_dat[[method]]

    # Rename columns
    for (old_name in names(rep_df)) {
      new_name <- ifelse(old_name %in% names(rename_mappings), paste0(method, rename_mappings[old_name]), old_name)
      names(rep_df)[names(rep_df) == old_name] <- new_name
    }

    # Add the dataframe to rep_dfs1
    rep_dfs1[[method]] <- rep_df
  }

  # Remove NULL entries
  rep_dfs1 <- rep_dfs1[!sapply(rep_dfs1, is.null)]

  # Subsetting those who are supposed to continue to Phase 2
  p2_ratio_eval = full[["true_P2_eval"]][[sim_lab]][["true_ratio_ind"]][["ratio_ind"]]
  index <- which(p2_ratio_eval == 0)

  # # Combine datasets
  # combined_data <- cbind(
  #   rep_dfs1[["czmk"]] %>%
  #     filter(subj.id %in% index) %>%
  #     dplyr::select(subj.id, czmk_OS, czmk_action, czmk_CIF) %>%
  #     group_by(czmk_action) %>%
  #     summarize(czmk_OS = mean(czmk_OS),
  #               czmk_CIF = mean(czmk_CIF)),
  #   rep_dfs1[["csk"]] %>%
  #     filter(subj.id %in% index) %>%
  #     dplyr::select(csk_OS, csk_action, csk_CIF)%>%
  #     group_by(csk_action) %>%
  #     summarize(csk_OS = mean(csk_OS),
  #               csk_CIF = mean(csk_CIF)),
  #   rep_dfs1[["zom"]] %>%
  #     filter(subj.id %in% index) %>%
  #     dplyr::select(zom_OS, zom_action, zom_CIF)%>%
  #     group_by(zom_action) %>%
  #     summarize(zom_OS = mean(zom_OS),
  #               zom_CIF = mean(zom_CIF))
  # )

  # Initialize an empty list to store the combined and summarized dataframes
  combined_data_list <- list()

  # Loop over each method
  for (method_ind in 1:length(methods)) {
    method = methods[method_ind]
    # print(method)

    if (!is.null(rep_dfs1[[method]])){
      # Subset the dataframe based on the index
      df <- rep_dfs1[[method]] %>%
        filter(subj.id %in% index)

      # Select columns and calculate means
      summary_df <- df %>%
        dplyr::select(subj.id, matches(paste0(method, "_"))) %>%
        group_by(!!sym(paste0(method, "_action"))) %>%
        summarize(
          !!sym(paste0(method, "_OS")) := mean(!!sym(paste0(method, "_OS"))),
          !!sym(paste0(method, "_CIF")) := mean(!!sym(paste0(method, "_CIF")))
        )

      if (method == "zom"){
        opposite = (1-df$zom_action[1])-1
        summary_df = rbind(c(opposite,0,0), summary_df)
      }

      if (method == "czmk"){
        summary_df = cbind(rep = sim,
                           tid = 1:2,
                           summary_df)
        method_labs = 1
      } else{
        method_labs = c(method_labs, method_ind)
      }

    } else{
      # print("null")
      summary_df = NA
    }


    # Store the summarized dataframe in the list
    combined_data_list[[method]] <- summary_df
  }

  # Combine the dataframes
  combined_data <- bind_cols(combined_data_list)

  combined_data_func <- function(var){
    combined_data_var = combined_data %>%
      dplyr::select(rep, tid, contains(var)) %>%
      pivot_longer(cols = contains(var),
                   names_to = c("method"),
                   values_to = var) %>%
      mutate(method = sub("_.*", "", method))
    return(combined_data_var)
  }

  combined_data_OS <- combined_data_func("OS")
  combined_data_CIF <- combined_data_func("CIF")
  combined_data_action <- combined_data_func("action") %>%
    mutate(action = factor(action))

  combined_data_long0 = left_join(combined_data_action, combined_data_OS,
                                  by = c("rep","tid", "method"))
  combined_data_long[[sim]] = left_join(combined_data_long0, combined_data_CIF, by = c("rep","tid", "method")) %>%
    as.data.frame() %>%
    mutate(method = factor(method,
                           levels = methods))
}

stacked_data = do.call(rbind, combined_data_long) %>%
  filter(!(method == "zom" & OS == 0 & CIF == 0)) %>%
  dplyr::select(-tid)
# Now stacked_data_filtered contains the filtered data

stacked_data = as.data.frame(stacked_data) %>%
  mutate(ncauses = fn$ncauses[i],
         # cause1prob = fn$cause1prob[i], # this only exists if endpoint = CR
         beta = fn$beta[i],
         prop = fn$prop[i],
         n    = fn$n[i],
         crit = fn$critS.no[i]) %>%
  mutate(
    method = factor(method,
                    levels = method.nm.simple,
                    labels = method.nm.abc),
    ncauses = factor(ncauses,
                     levels = 1:max(ncauses),
                     labels = c("1" = "2 causes")),
    setting = factor(beta,
                     levels = beta.levels,
                     labels = beta.labels),
    n = factor(n,
               levels = n.levels,
               labels = n.labels),
    design = factor(prop,
                    levels = design.levels,
                    labels = design.labels),
    crit.label = factor(crit,
                        levels = 1:2,
                        labels = c("1" = sprintf("Mean Truncated %s", crit_lab),
                                   "2" = sprintf("%s probability at t=%s, Curve(%s)", crit_lab, t0, t0)
                        )))

comb_plot = function(dataset,var){
  library(cowplot)
  if (var == "OS"){
    lab = "Overall Survival"
  } else{
    lab = "Cumulative Incidence Function"
  }
  ggplot(dataset, aes(x = method, y = !!as.name(var), color = action)) +
    geom_point() +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Add a boxplot without outliers
    labs(x = "Methods", y = "Probability", title = lab, color = "Treatment") +
    ylim(0,0.75) +
    scale_x_discrete(labels = paste0(method.nm.abc[method_labs], ": ", method.nm.formal[method_labs])) +
    theme_bw() +
    facet_grid(n ~ design + setting) +
    theme(legend.position = "bottom")
}
p1 = comb_plot(stacked_data, "OS")
p2 = comb_plot(stacked_data, "CIF")
p11 = p1 +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
p22 = p2 +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(y = NULL)
p.grid <- plot_grid(p11,p22,
                    align = "vh", axis = "tblr",
                    nrow = 1, ncol = 2,
                    common.legend = TRUE)
p.grid1 <- plot_grid(p.grid,get_legend(p1),
                     align = "vh",
                     axis = "tblr",
                     nrow = 2,
                     rel_heights = c(1, 0.1))

# Create an overarching title
overarching_title <- ggdraw() +
  draw_label(sprintf("FigX. Treatment Means across %s simulations for those whose diff in mean OS curve is within %s%%",
                     n.sim,
                     tol1),
             size = 10)
# Arrange the overarching title and combined plots vertically
final_plot <- plot_grid(overarching_title, p.grid1, ncol = 1, rel_heights = c(0.1, 0.9))
print(final_plot)
file.name.saved_p2 = file_naming(lab.date, "SimulationResults_P2Subset", crit.no)
save_plot(file.name.saved_p2, final_plot, base_height = 10, base_width = 20)

print("End of SM02.Plots_For_P2Subset_byTrt.R")
