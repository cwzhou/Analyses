# In this script, we calculate the TRUE OS curve for
# EVERYONE in the evaluation set comparing Trt -1 to Trt 1
# and identify which subjects have area under the curve
# within the pre-specified ratio
# NOTE: This is currently ONLY for AUC AKA MEANS ONLY

true1 = data.frame(id = 1:nrow(rep_czmk),
                   ratio = NA,
                   ratio_ind = NA)
true2 = list()
time_full = tp_surv
# below is copied from CR00.Simulation_Parameters.R
if (generate_failure_method == "simple_exp"){
  beta1.hazard0 = c(0, -1.6,-1.2,0.5)
  beta1.hazard1 = c(0, 0.3,-0.4,-0.2)
  beta2.hazard0 = c(0, 1.1,-1.3,0.3)
  beta2.hazard1 = c(0, -0.6,0.5,0.3)
  p = 0
} else{
  beta1.hazard0 = c(0,-1,-1.4)
  beta1.hazard1 = c(0,0.8,0.7)
  beta2.hazard0 = c(0,-0.2,1.2)
  beta2.hazard1 = c(0,-0.3,-2)
  p = cause1_prob$small.cause1prob$cause1prob
}

eps0_params = mean_tol1[1]

method.nm.abc =
  c("A", "B", "C")
method.nm.simple =
  c("czmk", "csk", "zom")
method.nm.formal =
  c("the proposed method", "Cho et al (2022)", "zero-order model")


os_simple_exp = function(cause1_prob, lam_cause1, lam_cause2, time_full){
  if (generate_failure_method == "simple_exp"){
    St = exp(-as.numeric((lam_cause1 + lam_cause2))*time_full)
  } else{
    term1 = (1-cause1_prob*(1-exp(-time_full)))^lam_cause1
    term2 = -(1-exp(-time_full*lam_cause2))*((1-cause1_prob)^lam_cause1)
    St = term1 + term2
  }
  return(St)
}

# this is for the eval set
for (person in 1:nrow(rep_czmk)){
  z1 = rep_czmk[person,]$Z1
  z2 = rep_czmk[person,]$Z2
  z3 = rep_czmk[person,]$Z3

  lamb_cause1_trt0 = as.numeric(exp(beta1.hazard0 %*% c(1,z1,z2,z3)))
  lamb_cause2_trt0 = as.numeric(exp(beta2.hazard0 %*% c(1,z1,z2,z3)))
  lamb_cause1_trt1 = as.numeric(exp(beta1.hazard1 %*% c(1,z1,z2,z3)))
  lamb_cause2_trt1 = as.numeric(exp(beta2.hazard1 %*% c(1,z1,z2,z3)))
  lambda_trt0 = c(lamb_cause1_trt0, lamb_cause2_trt0)
  lambda_trt1 = c(lamb_cause1_trt1, lamb_cause2_trt1)
  # taking the integration of OS curve from 0 to tau
  overs_t1 = integrate(os_simple_exp, lower = 0 , upper = tau,
                       cause1_prob = p,
                       lam_cause1 = lamb_cause1_trt1, lam_cause2 = lamb_cause2_trt1)$value
  overs_t0 = integrate(os_simple_exp, lower = 0 , upper = tau,
                       cause1_prob = p,
                       lam_cause1 = lamb_cause1_trt0, lam_cause2 = lamb_cause2_trt0)$value

  # taking the ratio of the means
  ratio = min(overs_t1/overs_t0, overs_t0/overs_t1)
  ratio_ind = ifelse(ratio > 1-eps0_params, 0, 1)
  true1["ratio"][person,] = ratio
  true1["ratio_ind"][person,] = ratio_ind

  # calculating the exact OS curve for a population with the same covariates as ID = person
  # for time_full points time_full
  os_trt1 = os_simple_exp(cause1_prob = p, lamb_cause1_trt1, lamb_cause2_trt1, time_full)
  os_trt0 = os_simple_exp(cause1_prob = p, lamb_cause1_trt0, lamb_cause2_trt0, time_full)
  true_OS_curve <- data.frame(x = time_full, y0 = os_trt0, y1 = os_trt1)

  plot = ggplot(true_OS_curve, aes(x = x)) +
    geom_line(aes(y = y0, color = "Trt-1"), linetype = "solid") +
    geom_line(aes(y = y1, color = "Trt1"), linetype = "dashed") +
    labs(x = "time_full", y = "St",
         title = sprintf("True OS Curve for a population the same as ID=%s", person),
         color = "Treatment") +
    scale_color_manual(values = c("Trt-1" = "purple", "Trt1" = "light blue"),
                       labels = c("Trt -1", "Trt 1")) +
    theme_minimal()

  true2[["beta1.hazard0"]] = beta1.hazard0
  true2[["beta1.hazard1"]] = beta1.hazard1
  true2[["beta2.hazard0"]] = beta2.hazard0
  true2[["beta2.hazard1"]] = beta2.hazard1
  true2[["eps0"]] = eps0_params
  true2[["covariates"]][person] = list(c(z1,z2,z3))
  true2[["true_OS_curve"]][person] = list(true_OS_curve)
  true2[["plot"]][person] = list(plot)
  true2[["lambda_trt0"]][person] = list(lambda_trt0)
  true2[["lambda_trt1"]][person] = list(lambda_trt1)
}

true2[["true_ratio_ind"]] = as.data.frame(true1)

names(rep_czmk)[names(rep_czmk) == "OS_eval"] <- "czmk_OS"
names(rep_csk)[names(rep_csk) == "OS_eval"] <- "csk_OS"
names(rep_zom)[names(rep_zom) == "OS_eval"] <- "zom_OS"
names(rep_czmk)[names(rep_czmk) == "CIF_eval"] <- "czmk_CIF"
names(rep_csk)[names(rep_csk) == "CIF_eval"] <- "csk_CIF"
names(rep_zom)[names(rep_zom) == "CIF_eval"] <- "zom_CIF"
names(rep_czmk)[names(rep_czmk) == "action"] <- "czmk_action"
names(rep_csk)[names(rep_csk) == "action"] <- "csk_action"
names(rep_zom)[names(rep_zom) == "action"] <- "zom_action"

# subsetting those who are supposed to continue to Phase 2
index = which(true1$ratio_ind == 0); length(index)
# Combine datasets
combined_data <- bind_cols(
  rep_czmk %>% filter(subj.id %in% index) %>% dplyr::select(subj.id,czmk_OS, czmk_action, czmk_CIF),
  rep_csk %>% filter(subj.id %in% index) %>% dplyr::select(csk_OS, csk_action, csk_CIF),
  rep_zom %>% filter(subj.id %in% index) %>% dplyr::select(zom_OS, zom_action, zom_CIF))

combined_data_func <- function(var){
  combined_data_var = combined_data %>%
    dplyr::select(subj.id, contains(var)) %>%
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

combined_data_long0 = left_join(combined_data_action, combined_data_OS, by = c("subj.id", "method"))
combined_data_long = left_join(combined_data_long0, combined_data_CIF, by = c("subj.id", "method")) %>%
  as.data.frame() %>%
  mutate(method = factor(method,
                         levels = c("czmk", "csk", "zom")))

rep_czmk %>% 
  filter(subj.id %in% index) %>%
  group_by(czmk_action) %>%
  summarize(n = n(),
            meanOS = mean(czmk_OS),
            meanCIF = mean(czmk_CIF))

rep_csk %>% 
  filter(subj.id %in% index) %>%
  group_by(csk_action) %>%
  summarize(n = n(),
            meanOS = mean(csk_OS),
            meanCIF = mean(csk_CIF))

calculate_statistics <- function(dataset, index_column, action_column, OS_column, CIF_column) {
  result <- dataset %>%
    filter({{index_column}} %in% index) %>%
    group_by({{action_column}}) %>%
    summarize(n = n(),
              meanOS = mean({{OS_column}}),
              meanCIF = mean({{CIF_column}}))
  return(result)
}

# Example usage with rep_czmk dataset
(result_czmk <- calculate_statistics(rep_czmk, subj.id, czmk_action, czmk_OS, czmk_CIF))
# Example usage with rep_csk dataset
(result_csk <- calculate_statistics(rep_csk, subj.id, csk_action, csk_OS, csk_CIF))
# Example usage with rep_zom dataset
(result_zom <- calculate_statistics(rep_zom, subj.id, zom_action, zom_OS, zom_CIF))
# Repeat for other datasets as needed




comb_plot = function(dataset,var){
  library(cowplot)
  if (var  == "OS"){
    lab = "Overall Survival"
  } else{
    lab = "Cumulative Incidence Function"
  }
  ggplot(dataset, aes(x = method, y = !!as.name(var))) +
    geom_point(aes(color = action)) +
    geom_boxplot(fill = "#E0BBE4", color = "#673AB7",alpha = 0.2, outlier.shape = NA) +  # Add a boxplot without outliers
    labs(x = "Methods", y = "Probability", title = lab, color = "Treatment") +
    # ylim(0.6,1.2) + # this needs to change depending on the generate_failure_time setting
    scale_x_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}
p1 = comb_plot(combined_data_long, "OS")
p2 = comb_plot(combined_data_long, "CIF")
p.1 = p1 +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
p.2 = p2 +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(y = NULL)
p.grid <- plot_grid(p.1,p.2,
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
  draw_label(sprintf("For those whose OS curve is within %s%% (N = %s)",
                     eps0_params, length(index)),
             size = 10)
# Arrange the overarching title and combined plots vertically
final_plot <- plot_grid(overarching_title, p.grid1, ncol = 1, rel_heights = c(0.1, 0.9))
print(final_plot)


combined_data_long %>% group_by(method) %>%
  summarise(n = n(),mean_os = mean(OS), mean_cif = mean(CIF))
combined_data_long %>% group_by(method, action) %>%
  summarise(n = n(),mean_os = mean(OS), mean_cif = mean(CIF))
