# sbatch -p general -N 1 --mem=50GB -n 2 -t 10-07:00:00 --mail-type=end --mail-user=cwzhou@email.unc.edu --wrap="Rscript CR01.Simulation_Run.R > /work/users/c/w/cwzhou/Proj1/CRoutput_20250126_test.txt"
# For local: install.packages('~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/itrSurv_0.1.0.tar.gz', repos = NULL, type = 'source')
# For cluster: install.packages('/nas/longleaf/home/cwzhou/Dissertation/itrSurv/itrSurv_0.1.0.tar.gz', repos = NULL, type = 'source')
# bash CR_S2value.sh

## This is the script used to run the simulations ##
## Edit parameters in CR00.Simulation_Parameters.R script
## Plots in CR02.Simulation_Summary.R script
# setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
local = 1 # change to 0 for cluster
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
} else{
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper1_CR")
}

# call on parameters which calls on library
source("CR00.Simulation_Parameters.R")

# call on functions then runs simulations
if (parallel == 1){
  message("WE ARE RUNNING PARALLEL JOBS.")
  source("CR00.Simulation_Body.R")
  print("now going into rbind")
  # Combine the results into a single dataframe
  final_results0 <- do.call(rbind, results_list) # years
  # print(is.data.frame(final_results0))
  final_results0 = as.data.frame(final_results0)
  # print(is.data.frame(final_results0))
  # print(tail(final_results0))
  result_sub0 = final_results0 %>%
    dplyr::select(obs_survival, czmk_survival,
                  csk_survival, pmcr_survival,
                  aipwe_survival, zom_survival,
                  obs_endpoint, czmk_endpoint,
                  csk_endpoint, pmcr_endpoint,
                  aipwe_endpoint, zom_endpoint)
  means0 = apply(result_sub0, 2, mean, na.rm=TRUE)
  sds0 = apply(result_sub0, 2, sd, na.rm=TRUE)
  mean_sd0 = cbind(means0,sds0)
  
  final_results <- final_results0 %>%
    #   # convert to days
    mutate(across(ends_with("_survival") | ends_with("_endpoint"), ~ . * 365.25))
  print("final_results")
  print(is.data.frame(final_results))
  print(tail(final_results,10))
} else{
  message("Not running parallel")
  #not parallel
  source("CR00.Simulation_Body_noparallel.R")
  final_results = result
}


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

  # message('saving true2')
  # saveRDS(list(settings = setting,
  #              true2_full_P2_eval = true2,
  #              ind_stat = trt_result),
  #         filename.true2)
}

# source("plotting value functions.R")
# plot_values(combined_data, OS_eval)
# plot_values(combined_data, CIF_eval)

if (local == 1){
  # beep()
  # data.df%>%group_by(action) %>% summarize(mean = mean(event.time))
  # data.df%>%group_by(status,action) %>% summarize(mean = mean(event.time))
  print(final_results)
  print(round(mean_sd,3))
  # View(final_results)

  # if (!skip.csk){
  #   print(rep_csk %>%
  #           group_by(action) %>%
  #           summarise(n=n(),
  #                     os = mean(OS_eval),
  #                     cif = mean(CIF_eval)))
  # }
  # 
  # if (!skip.czmk){
  #   print(rep_czmk %>%
  #           group_by(action) %>%
  #           summarise(n=n(),
  #                     os = mean(OS_eval),
  #                     cif = mean(CIF_eval)))
  # }
}

message("End of CR01.Simulation_Run.R -- proceed to CR02.Simulation_Summary.R for creating plots.")
