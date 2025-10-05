# Record start time
start_time <- Sys.time()
# bash CR_S2value.sh

## This is the script used to run the simulations ##
## Edit parameters in CR00.Simulation_Parameters.R script
## Plots in CR02.Simulation_Summary.R script

# single machine/multi-node (using CR_S2value.sh)
# call on parameters which calls on library
source("CR00.Simulation_Parameters.R")
print(sim_deets)
source("CR00.Simulation_Body_noparallel_blindedreview.R")
final_results = result

result_sub = final_results %>%
  dplyr::select(obs_survival, czmk_survival,
                csk_survival, pmcr_survival,
                aipwe_survival, zom_survival,
                obs_endpoint, czmk_endpoint,
                csk_endpoint, pmcr_endpoint,
                aipwe_endpoint, zom_endpoint,
                obs_trt1, czmk_trt1,
                csk_trt1, pmcr_trt1,
                aipwe_trt1, zom_trt1,
                czmk_n_phase2,
                zom_n_phase2,
                training_percent.censor)
means = apply(result_sub, 2, mean, na.rm=TRUE)
sds = apply(result_sub, 2, sd, na.rm=TRUE)
mean_sd = cbind(means,sds)
rounded_mean_sd = round(mean_sd,3)

long_res = list(statistics = final_results,
                settings = setting,
                # true_P2_eval = true3,
                # p2.df = p2.df,
                mean_sd = mean_sd,
                rounded_mean_sd = rounded_mean_sd)

#not saving result_ext right now (statistics = result not statistics = result_ext)
if (savingrds == TRUE){
  message('saving rds')
  saveRDS(long_res,
          filename)
}

if (local == 1){
  print(final_results)
  print(round(mean_sd,3))
}

# Record end time
end_time <- Sys.time()
# Calculate total duration
total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
hours <- floor(total_secs / 3600)
minutes <- floor((total_secs %% 3600) / 60)
seconds <- round(total_secs %% 60)
# Print results in readable format
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

print(sim_deets)

cat("Total time:", hours, "hours", minutes, "minutes", seconds, "seconds\n")
message("End of CR01.Simulation_Run.R -- proceed to CR02.Simulation_Summary.R for creating plots.")

