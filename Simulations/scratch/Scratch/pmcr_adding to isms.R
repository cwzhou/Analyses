# estimation
cat ("4. pmcr \n")
if (!skip.pmcr) {
  cat ("  4. pmcr - Policy estimation \n")
  arg.pmcr2 = list(Time="obs_time",
                   Event="status",
                   formula=Trt ~ Z1 + Z2, # Trt has to be 0/1 for pmcr
                   data=data.df.pmcr, # b/c trt is 0/1
                   rgenoud=FALSE, 
                   Restrict=FALSE, # required for unrestricted
                   propscore="logistic", # change as appropriate
                   t0=1) # t0 need to be tuned
  # Unrestricted regime
  set.seed(train_seed + 3)
  temp.unrestr.fit <- do.call(PMCR, arg.pmcr2)
  #range of alpha # alpha needs to be tuned
  alps<-c(temp.unrestr.fit$Fbeta2[2],
          seq(round(temp.unrestr.fit$Fbeta2[2],2)+0.01,
              round(temp.unrestr.fit$Fbeta1[2],2)+0.03,
              0.01)) 
  # NEED TO DO: CODE TUNING FOR ALPHA; update t0 above to be appropriate as well
  #Restricted regime parameters
  arg.pmcr3 = arg.pmcr2; 
  arg.pmcr3$Restrict = TRUE
  alp<-0.1; arg.pmcr3$alp = alp # alp needs to be tuned based on alps
  M<-1e+5; arg.pmcr3$M = M # pick a large M
  #Restricted regime
  set.seed(train_seed + 3)
  temp.restr.fit <- do.call(PMCR, arg.pmcr3)
  optimal.pmcr <- temp.restr.fit$beta3 #optbetas for the linear decision rule

  pmcr.error <- class(optimal.pmcr)[1] == "try-error"
  arg.pmcr$policy <- if (!pmcr.error) optimal.pmcr
  policy_pmcr <<- arg.pmcr$policy
  rm(optimal.pmcr); gc()
  
  cat ("  \n 4. pmcr - Evaluation \n")
  arg_pmcr <<- arg.pmcr
  set.seed(train_seed + 10)
  if (!pmcr.error) {
    pmcr.data.rep <- do.call(gdata_CR, arg.pmcr)
    rep_pmcr <<- pmcr.data.rep
    result[sim, "pmcr_survival"] = overall_survival_val.fn(pmcr.data.rep)
    # overall_survival_mean <- mean(subset(pmcr.data.rep$event.time, pmcr.data.rep$status > 0),
    #                               na.rm = TRUE)
    # endpoint_mean <- mean(subset(pmcr.data.rep$event.time, pmcr.data.rep$status == priority_cause),
    #                       na.rm = TRUE)
    result[sim, "pmcr_endpoint"] = endpoint_val.fn(pmcr.data.rep, priority_cause)
    result[sim, "pmcr_percent.censor"] <- mean(pmcr.data.rep$status==0, na.rm = TRUE)
    result[sim, paste0("pmcr_cause.", 1:n.causes)] <-
      sapply(1:n.causes, function(s) mean(pmcr.data.rep$status == s))
    # if (!pmcr.error) result[sim, "pmcr_phase1"] <- val.fn_phase1(pmcr.data.rep$event.time)
    # if (!pmcr.error) result[sim, "pmcr_phase2"] <- val.fn_phase2(pmcr.data.rep$event.time)
  }
  # within pmcr method, mean survival time from those who have action -1 and mean survival time for those who have action 1
  A_mean = pmcr.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
  B_mean = pmcr.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
  A_mean_cr = pmcr.data.rep %>% filter(action == -1) %>% endpoint_val.fn(priority_cause)
  B_mean_cr = pmcr.data.rep %>% filter(action == 1) %>% endpoint_val.fn(priority_cause)
  trt_result[sim_count+1, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("pmcr", A_mean, B_mean, A_mean_cr, B_mean_cr)
  rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
  
  result[sim, "time.pmcr"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
  arg.pmcr$policy <- NULL; gc()
  rm(pmcr.data.rep); gc()
  
  # temp.unrestr.fit<-PMCR(Time="obs_time",
  #                        Event="status",
  #                        formula=Trt ~ Z1 + Z2,
  #                        data=data.df.pmcr,
  #                        rgenoud=FALSE,
  #                        Restrict=FALSE,
  #                        propscore="logistic",
  #                        t0=1)
  # temp.restr.fit<-PMCR(Time="time",
  #                      Event="status",
  #                      formula=modelPr1,
  #                      data=train,
  #                      rgenoud=FALSE,
  #                      Restrict=TRUE,
  #                      propscore="logistic",
  #                      t0=t0,
  #                      alp=alp,
  #                      M=M)
}