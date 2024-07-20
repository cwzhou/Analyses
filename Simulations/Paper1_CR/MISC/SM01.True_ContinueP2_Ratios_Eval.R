
# right now: this only works for crit = mean or area (NOT prob)
# right now: this only works for simple_exp setting
true1 = data.frame(id = 1:n.eval,
                   ratio = NA,
                   ratio_ind = NA)
time_full = tp_surv
# # below is copied from CR00.Simulation_Parameters.R
# beta1.hazard0 = c(0, -1.6,-1.2,0.5)
# beta1.hazard1 = c(0, 0.3,-0.4,-0.2)
# beta2.hazard0 = c(0, 1.1,-1.3,0.3)
# beta2.hazard1 = c(0, -0.6,0.5,0.3)
eps0_params = mean_tol1[1]
p = setting[["cause1prob"]]

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

true2[["setting"]] = setting
# true2[["beta1.hazard0"]] = beta1.hazard0
# true2[["beta1.hazard1"]] = beta1.hazard1
# true2[["beta2.hazard0"]] = beta2.hazard0
# true2[["beta2.hazard1"]] = beta2.hazard1
true2[["eps0"]] = eps0_params
# true2[[sim]][["tp_surv"]] = list(tp_surv)

# this is for the eval set
for (person in 1:n.eval){

  # Extract the relevant row for the person
  person_row <- rep_czmk[person, ]
  # Extract the 'z' columns
  z_columns <- person_row[grep("^Z", names(person_row))]
  # List of beta hazards
  beta_hazards <- list(beta1.hazard0, beta2.hazard0, beta1.hazard1, beta2.hazard1)
  results <- list()
  # Loop through beta hazards and calculate corresponding lambdas
  for (i in seq_along(beta_hazards)) {
    results[[i]] <- as.numeric(exp(beta_hazards[[i]] %*% c(1, as.numeric(z_columns))))
  }

  # Assign results to variables
  lamb_cause1_trt0 <- results[[1]]
  lamb_cause2_trt0 <- results[[2]]
  lamb_cause1_trt1 <- results[[3]]
  lamb_cause2_trt1 <- results[[4]]
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

  # below takes too much memory so only do this if NEEDED
  # # calculating the exact OS curve for a population with the same covariates as ID = person
  # # for time_full points time_full
  os_trt1 = os_simple_exp(p, lamb_cause1_trt1, lamb_cause2_trt1, time_full)
  os_trt0 = os_simple_exp(p, lamb_cause1_trt0, lamb_cause2_trt0, time_full)
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

  sim_lab = sprintf("sim%s", sim)
  # uncomment below if plots are run
  true2[[sim_lab]][["covariates"]][person] = list(z_columns)
  # true2[[sim_lab]][["true_OS_curve"]][person] = list(true_OS_curve)
  # true2[[sim_lab]][["plot"]][person] = list(plot)
  true2[[sim_lab]][["lambda_trt0"]][person] = list(lambda_trt0)
  true2[[sim_lab]][["lambda_trt1"]][person] = list(lambda_trt1)
  true2[[sim_lab]][["true_ratio"]][person] = list(ratio)
  true2[[sim_lab]][["true_ratio_ind"]][person] = list(ratio_ind)
}

rep_dfs <- list()
rep_dfs <- list(
  "czmk" = if (!get("skip.czmk", envir = .GlobalEnv)) rep_czmk else NULL,
  "csk" = if (!get("skip.csk", envir = .GlobalEnv)) rep_csk else NULL,
  "pmcr" = if (!get("skip.pmcr", envir = .GlobalEnv)) rep_pmcr else NULL,
  "zom" = if (!get("skip.zom", envir = .GlobalEnv)) rep_zom else NULL,
  "obs" = if (!get("skip.obs", envir = .GlobalEnv)) rep_obs else NULL
)

true3[[sim_lab]][["true_ratio_ind"]] = as.data.frame(true1)
true3[[sim_lab]][["eval_datasets"]] = rep_dfs

message("End of SM01.True_ContinueP2_Ratio_Eval.R")
