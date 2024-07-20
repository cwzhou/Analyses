library(dplyr);library(ggplot2);library(cowplot); library(tidyr)

testing_out = 1

crit.no = 1
critS.no = 1
critE.no = 1
endpoint = "CR"

lab.date = "2023-11-27" # change this for the date of RDS data you want
files <- list.files(path = sprintf("%s/output", endpoint), pattern = paste0(lab.date, ".*\\.rds"), full.names = TRUE)
if (testing_out == 1){
  files = files[grepl("_tmp", files) == FALSE]
  print(files)
}
method.nm.abc =
  c("A", "B", "C") #, "D", "E", "F")
method.nm.simple =
  c("czmk", "zom", "observed")
# "czmk_endpoint", "zom_endpoint", "observed_endpoint") #"csk", "gkRF", "gkLM", "dw",
method.nm.formal =
  c("the proposed method", "zero-order model", "observed policy")#,
# "the proposed method - ep", "zero-order model - ep", "observed policy - ep")
# "Cho et al. (2022)", "Goldberg & Kosorok (2012), RF",
# "Goldberg & Kosorok (2012), linear", "Simoneau et al. (2019)",

filename = function(lab.date, crit.no, Phase){
  paste0(endpoint,"/figure/C22_", Phase, "Summary_", gsub("-", "", lab.date), "_crit", crit.no, ".eps")
}

if (endpoint == "CR"){
  fn <-
    expand.grid(ncauses = 1, cause1prob = 1:2,
                beta = 1, prop = 1:2, n = 1:2, critS.no = 1:2, critE.no = 1:2) %>%
    mutate(nm = paste0(beta, "-", prop, "-", n),
           fn = paste0(endpoint,"/output/simResult_", lab.date,
                       "_nCauses", ncauses, "_cause1prob", cause1prob,
                       "_beta", beta, "_prop", prop,
                       "_n", n, "_critS", critS.no, "_critE", critE.no, ".rds"))
} else{
  fn <-
    expand.grid(ncauses = 1, beta = 1:6, prop = 1:2, n = 1:2, critS.no = 1:2, critE.no = 1:2) %>%
    mutate(nm = paste0(beta, "-", prop, "-", n),
           fn = paste0(endpoint,"/output/simResult_", lab.date,
                       "_nCauses", ncauses,
                       "_beta", beta, "_prop", prop,
                       "_n", n, "_critS", critS.no, "_critE", critE.no, ".rds"))
}

for (Phase in 1:2){
  if (Phase == 1){
    Phase_lab = "survival"
  } else{
    Phase_lab = "endpoint"
  }
  method.nm.simple1 = paste0(method.nm.simple, "_", Phase_lab)
  result.comb <-
    lapply(1:dim(fn)[1], function(i) {
      fn.i = fn$fn[i]
      if (file.exists(fn.i)) {
        print(sprintf("file %s exists", fn.i))
        full <- readRDS(fn.i)
        a <- full$statistics
      } else if (file.exists(gsub("\\.rds", "_tmp.rds", fn.i))) {
        print(sprintf("_TMP.RDS exists: %s", fn.i))
        full <- NULL
        a <- readRDS(gsub("\\.rds", "_tmp.rds", fn.i))
      } else {
        full <- NULL
        a <- NULL
      }
      if (is.null(full)){
        NULL
      } else{
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
        }
    }
    )}