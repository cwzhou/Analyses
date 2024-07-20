# arg3 beta
if (endpoint == "CR"){
  
  #   print("simple exp coeff")
  # betas <- list(
  #   beta1 = list(
  #     beta1.hazard0 = c(0,-0.2,-0.5),
  #     beta1.hazard1 = c(0,0.1,0.4),
  #     beta2.hazard0 = c(0,-0.4,-0.7), #c(0,-0.1,-0.2),
  #     beta2.hazard1 = c(0,-1.1,-1.5))
  # )

  # zom breaks down here
  betas <- list(
    beta1 = list(
      beta1.hazard0 = c(0,-0.4,-0.5,0.3),
      beta1.hazard1 = c(0,0.1,0.4,-0.1),
      beta2.hazard0 = c(0,-0.4,-0.7,0.6), #c(0,-0.1,-0.2),
      beta2.hazard1 = c(0,-1.1,-1.5,-1))
  )

  alpha1 <- 0 #log(lambda0)
  alpha2 <- 0 #log(lambda0)
  
  ncov.list <- lapply(betas, function(x) length(x$beta1.hazard1) -1 ) # removing intercept
} else{
  stop("betas not coded up yet for endpoint: ", endpoint)
}

# setting1 n.boot 50 / n 300
# this setting is for the first of each thing, since arg = all 1s
# args 1- 7: "endpoint", "ncauses", "beta", "propensity",
# "size", "crit_surv", "crit_endpoint"
if (endpoint == "CR"){
  setting = c(all_methods = list(all_methods),
              n.methods = n.methods,
              arg = list(arg), default, ncauses[[arg2]], cause1_prob[[arg8]],
              betas[[arg3]], ncov = ncov.list[[arg3]],
              propensity[[arg4]], size[[arg5]],
              crit_surv[[arg6]], crit_endpoint[[arg7]])
} else{
  setting = c(all_methods = all_methods,
              n.methods = n.methods,
              arg = list(arg), default, ncauses[[arg2]],
              betas[[arg3]], ncov = ncov.list[[arg3]],
              propensity[[arg4]], size[[arg5]],
              crit_surv[[arg6]], crit_endpoint[[arg7]])
}


### 2. put a selected setting into the global environment
cat("setting (endpoint, ncauses, beta, propensity,
    size, crit_phase1, crit_phase2, cause1_prob) ", arg, "\n")
list2env(setting, envir = globalenv())
beta.propensity <- beta.propensity(ncov)
print('test3')
if (endpoint == "CR"){
  if (testing != 1){
    filename = paste0(endpoint,"/output/",generate_failure_method,"/simResult_", arg.date,
                      "_nCauses", arg2, "_cause1prob", arg8,
                      "_beta", arg3, "_prop", arg4,
                      "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")
  } else{
    filename = paste0(endpoint,"/output/",generate_failure_method,"/testing_simResult_", arg.date,
                      "_nCauses", arg2, "_cause1prob", arg8,
                      "_beta", arg3, "_prop", arg4,
                      "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")
  }
  
} else{
  filename = paste0(endpoint,"/output/",generate_failure_method,"/simResult_", arg.date,
                    "_nCauses", arg2,
                    "_beta", arg3, "_prop", arg4,
                    "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")
}
print(filename)

# complete the scores
predPropensityFn <- function(covariate) {
  # message("predPropensityFn")
  plogis(cbind(1, covariate) %*% beta.propensity)
}
if (endpoint == "CR"){
  # alpha1/2 correspond to baseline hazard (to make failure times higher)
  # alpha1=alpha2=0 for now; if i change that then need to include in "settings"
  predHazardFn <- function(action, covariate, cause) {

    if (cause == 1){
      message("cause1")
      ifelse(action == 1,
             alpha1 + cbind(1, covariate) %*% beta1.hazard1, #cause1 treatment1
             alpha1 + cbind(1, covariate) %*% beta1.hazard0) #cause1 treatment 0
      # mat_beta = c(beta11,beta12,beta13)
      # # print(dim(mat_cov))
      # # print(mat_beta)
      # # print(dim(mat_beta))
      # alpha1 + mat_cov %*% mat_beta
    } else if (cause == 2){
      message("cause2")
      ifelse(action == 1,
             alpha2 + cbind(1, covariate) %*% beta2.hazard1, #cause2 treatment1
             alpha2 + cbind(1, covariate) %*% beta2.hazard0) #cause2 treatment0
      # mat_beta = c(beta21,beta22,beta23)
      # alpha2 + mat_cov %*% mat_beta
    } else if (cause == 3){
      message("cause3")
      ifelse(action == 1,
             alpha2 + cbind(1, covariate) %*% beta3.hazard1, #cause3 treatment1
             alpha2 + cbind(1, covariate) %*% beta3.hazard0) #cause3 treatment0
      # mat_beta = c(beta21,beta22,beta23)
      # alpha2 + mat_cov %*% mat_beta
    }
    else{
      stop("CURRENTLY ONLY CODED FOR 2 (or 3?) CAUSES")
    }
  }
} else{
  predHazardFn <- function(action, covariate, cause = NULL) {
    ifelse(action == 1,
           cbind(covariate) %*% beta.hazard1,
           cbind(covariate) %*% beta.hazard0)
  }
}
