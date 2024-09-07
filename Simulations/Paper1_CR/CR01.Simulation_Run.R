## This is the script used to run the simulations ##
## Edit parameters in CR00.Simulation_Parameters.R script
## Plots in CR02.Simulation_Summary.R script
# setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
<<<<<<< HEAD
local = 1
=======
local = 0
>>>>>>> 6beef6b5be9d9621ad573e2380672c9f91880dea
parallel = 0
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
} else{
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper1_CR")
}

# call on parameters which calls on library
source("CR00.Simulation_Parameters.R")

# call on functions then runs simulations
if (parallel == 1){
  source("CR00.Simulation_Body.R")
  print("now going into rbind")
  # Combine the results into a single dataframe
  final_results0 <- do.call(rbind, results_list) # years
  print(is.data.frame(final_results0))
  final_results0 = as.data.frame(final_results0)
  print(is.data.frame(final_results0))
  print(tail(final_results0))
  final_results <- final_results0 %>%
    #   # convert to days
    mutate(across(ends_with("_survival") | ends_with("_endpoint"), ~ . * 365.25))
  print("final_results")
  print(is.data.frame(final_results))
  print(tail(final_results,10))
} else{
  source("CR00.Simulation_Body_noparallel.R")
  final_results = result
}


result_sub = final_results %>%
  dplyr::select(obs_survival, czmk_survival,
                csk_survival, pmcr_survival,
                aipwe_survival, zom_survival,
                obs_endpoint, czmk_endpoint,
                csk_endpoint, pmcr_endpoint,
                aipwe_endpoint, zom_endpoint)


means = apply(result_sub, 2, mean, na.rm=TRUE)
sds = apply(result_sub, 2, sd, na.rm=TRUE)
mean_sd = cbind(means,sds)

long_res = list(statistics = final_results,
                settings = setting,
                # true_P2_eval = true3,
                # p2.df = p2.df,
                mean_sd = mean_sd)

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
# write.csv(final_results, "/nas/longleaf/home/js9gt/survrf/Outputs/1STRATA_28Aug_10stage__500pt_lowcens_V2", row.names=FALSE)


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

  print(rep_csk %>%
          group_by(action) %>%
          summarise(n=n(),
                    os = mean(OS_eval),
                    cif = mean(CIF_eval)))

  print(rep_czmk %>%
          group_by(action) %>%
          summarise(n=n(),
                    os = mean(OS_eval),
                    cif = mean(CIF_eval)))
}

message("End of CR01.Simulation_Run.R -- proceed to CR02.Simulation_Summary.R for creating plots.")
