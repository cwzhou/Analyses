# full = readRDS(".CR/output/simResult_2023-11-27_nCauses1_cause1prob1_beta1_prop1_n1_critS2_critE2.rds")
library(xtable)
mean_sd_df = as.data.frame(full$mean_sd)

# Create and print tables for each method type
means_table = list()
for (type in c("survival", "endpoint")) {
  if (type == "endpoint"){
    type1 = "CIF"
    label = full$settings$value_phase2
  } else{
    type1 = type
    label = full$settings$value_phase1
  }
  # cat(paste0("Mean ", toupper(type1), " for ", full$settings$n.sim, " simulations\n"))
  cat(paste0("Mean and SD of ", tools::toTitleCase(sub("mean.*", "mean", label)), "s across ", full$settings$n.sim, " simulations\n"))
  cat(paste0("\nSimulation Setting:", 
             "\nNumber of Subjects in Training Set: ", full$settings$n,
             "\nNumber of Subjects in Testing Set: ", full$settings$n.eval,
             "\nNumber of Causes: ", full$settings$M,
             "\nCause1prob: ", full$settings$cause1prob,
             "\nNumber of Covariates: ", full$settings$ncov,
             "\nCovariate Coefficient for Cause 1 Treatment 0: beta1.hazard0: ", toString(full$settings$beta1.hazard0),
             "\nCovariate Coefficient for Cause 1 Treatment 1: beta1.hazard1: ", toString(full$settings$beta1.hazard1),
             "\nCovariate Coefficient for Cause 2 Treatment 0: beta2.hazard0: ", toString(full$settings$beta2.hazard0),
             "\nCovariate Coefficient for Cause 2 Treatment 1: beta2.hazard1: ", toString(full$settings$beta2.hazard1),
             "\nFirst Endpoint (Survival): ", full$settings$value_phase1,
             "\nSecond Endpoint (", full$settings$endpoint, "): ", full$settings$value_phase2,"\n","\n"
             ))
  
  # Filter data
  filtered_data <- mean_sd_df[grepl(paste0("_", type, "$"), rownames(mean_sd_df)), ]
  
  # Remove suffix from row names
  row_names <- sub(paste0("_", type, "$"), "", rownames(filtered_data))
  rownames(filtered_data) <- row_names
  
  # Print the table
  print(filtered_data)
  means_table[[type]] = filtered_data
  
  cat("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
}

surv_means_table_latex = xtable(means_table$survival,
                          caption = sprintf("Mean and Sd of Survival Endpoint for all Simulations (n.sim = %s)", full$settings$n.sim))
endpoint_means_table_latex = xtable(means_table$endpoint,
                               caption = "Mean and Sd of Survival Endpoint for all Simulations in this Setting")
# print(surv_means_table_latex, include.rownames = TRUE)
# print(endpoint_means_table_latex, include.rownames = TRUE)


# # next table
# ind_dat = full$statistics %>% dplyr::select(CZMKsurv_ind_vsZOM, CZMKcif_ind_vsZOM,
#                                            CZMKsurv_ind_vsOBS, CZMKcif_ind_vsOBS,
#                                            ZOMsurv_ind_vsOBS, ZOMcif_ind_vsOBS
#                                            )
# tab2 = as.data.frame(apply(ind_dat, 2, mean))
# colnames(tab2) <- c("mean")
# # tab2 is the proportion of simulations where first method mean value (of subjects) is within 1-eps0 of second method value (survival) or below second method (endpoint = CR)
# # mean of simulations (each sim value is mean over all subjects)
# # mean(full$statistics$CZMKsurv_ind_vsZOM)
