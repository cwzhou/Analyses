cat("\n******************************\n")
# estimation
cat ("4. aipwe for Simulation",sim, ":",generate_failure_method,"\n")
if (!skip.aipwe) {
  cat ("  4. aipwe - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
  # models_aipwe <- paste0("Trt ~ ",
                        # paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
  arg.aipwe2 = list(data=data.df.aipwe, # b/c trt is 0/1
                   pp.v = ncov,
                   tau1 = t0_aipwe, # t0 needs to be tuned
                   tune = c(0.001,0.01,0.5,1,
                            seq(0.1,400,
                                length.out=16))
                   ) 
  set.seed(train_seed + 4)
  aipwe.fit1 <- do.call(aipwe.fit, arg.aipwe2)
  # aipwe.fit(data_list = data.df.aipwe,
  #           pp.v = ncov,
  #           tau1 = 0.1,
  #           tune = tune)
  optimal.aipwe <- aipwe.fit1 #betas
  
  aipwe.error <- class(optimal.aipwe)[1] == "try-error"
  arg.aipwe$policy <- if (!aipwe.error) optimal.aipwe
  policy_aipwe <<- arg.aipwe$policy
  rm(optimal.aipwe); gc()
  
  cat ("  \n 4. aipwe - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
  arg_aipwe <<- arg.aipwe
  set.seed(train_seed + 10)
  if (!aipwe.error) {
    aipwe.data.rep <- do.call(gdata_CR, arg.aipwe)
    rep_aipwe <<- aipwe.data.rep
    result[sim, "aipwe_survival"] = overall_survival_val.fn(aipwe.data.rep)
    result[sim, "aipwe_endpoint"] = endpoint_val.fn(aipwe.data.rep)
  }
  # within aipwe method, mean survival time from those who have action -1 and mean survival time for those who have action 1
  A_mean = aipwe.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
  B_mean = aipwe.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
  A_mean_cr = aipwe.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
  B_mean_cr = aipwe.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
  trt_result[sim_count+3, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("aipwe", A_mean, B_mean, A_mean_cr, B_mean_cr)
  rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
  
  result[sim, "time.aipwe"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
  arg.aipwe$policy <- NULL; gc()
  rm(aipwe.data.rep); gc()
}