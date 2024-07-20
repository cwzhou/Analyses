    N=10;
    ncov = 2;
    mass_p = 0.3;
    predHazardFn; predPropensityFn;
    ztype=0; zparam=0.5;
    ctype=1;cparam=2;censor_min=0;censor_max=NULL;
    num_A=2;tau=5;
    seed1=2023
    n.phases = 2;
    tau = 10;          # structural parameters
    Sig = diag(ncov) + 0.2 - diag(ncov) * 0.2;      # covariate structure
    at.risk = rep(1, N);
    predHazardFn;      # list of predictor functions
    predPropensityFn;  # optimal rule (if !is.null; propensity scores are ignored.) for value calculation
    # policy = NULL     # policy is an object of rsf.obj list.


    if (ztype == 0){z <- matrix(rnorm(N * ncov), nrow = N, ncol = ncov)} #print("Covariates follow N(0,1)");
    if (ztype == 1){z <- matrix(rbinom(N*ncov, 1, zparam),nrow=N,ncol=ncov)}
    if (ztype == 2){z <- matrix(runif(N*ncov),nrow=N,ncol=ncov)}
    if (is.null(colnames(z))) colnames(z) = paste0("Z", 1:ncov)
    covariate <- z
    # generating censoring time
    if (ctype == 0){cc <- rexp(N,cparam)}
    if (ctype == 1){
    if(is.null(censor_max)){censor_max = tau}
    if(is.null(censor_min)){censor_min = 0}
    censor_time <- cc <- runif(N,min=censor_min,max=censor_max)}

    message("multiphasedynamics")
    library(MASS)
    nphases = n.phases
    require(dplyr)
    if (length(at.risk) == 1) at.risk = rep(at.risk, N)
    if (length(at.risk) != N)  stop ("length of at.risk should match.")
    if (!is.null(covariate)) if (dim(covariate)[1] != N) stop ("length of at.risk should match.")

    # initial state
    if (is.null(covariate))     covariate = mvrnorm(N, rep(0, ncov), Sig)
    if (is.null(colnames(covariate))) colnames(covariate) = paste0("Z", 1:ncov)
    if (is.null(censor_time)){
      print("censor time was not generated earlier, randomly generating from unif(0,5) - change if needed")
      censor_time = runif(N,min=0,max=5)
    }

    # skeleton
    output <-
      array(NA, dim = c(N, 2 + 2 + 2 + ncov),
            dimnames = list(1:N, c("subj.id", "rep.id", "event.time", "status",
                                   "action", "at.risk", colnames(covariate))))
    output[, "subj.id"] <- 1:N
    output[, "rep.id"] <- 0           # rep.id is reserved for later use (repeated random draws)
    action = rep(NA, N)           # initialize action vector

    output[, "at.risk"]           <- at.risk
    output[, colnames(covariate)] <- covariate
    print(output)

    ## action                (t=0)
    if (!is.null(policy)) {
      if (attr(policy, "class") %in% c("ITRSurv")) { #### new package!
        print("policy is from ITRSurv estimation")
        x = output[as.logical(at.risk),] %>% output2observable()
        x$Trt = x$A
        tmp_act = matrix(nrow = N, ncol = nphases+1+1)
        colnames(tmp_act) = c("Stop", "P1", "P2", "Final")
        # below, we want D.1 b/c priority cuase?
        for (phase in 1:nphases){
          x[, c("obs_time", paste0("D.", phase-1), "Trt", "A")] = 0  # dummy values in order for the package not to drop the NA rows.
          x = get_all_vars(policy@phaseResults[[phase]]@model, x)
          args <- list(Phase = phase, policy, newdata = x)
          docalling <<- do.call(predict, args)
          if (phase == 1){
            stop = docalling$optimal@Ratio_Stopping_Ind #stop=1 for P1 and stop=0 for P2
            print(stop)
            tmp_act[,1] = stop
          }
          tmp_act[,phase+1] <- docalling$optimal@optimalTx
          print(tmp_act)
        }
        tmp_act[, 4] <- ifelse(tmp_act[, "Stop"] == 1,
                               tmp_act[, "P1"],
                               tmp_act[, "P2"])
        print(tmp_act)
        output[,"action"] = tmp_act[,4]
        action[at.risk != 0] <- do.call(predict, args)$optimal@optimalTx # action taken at phase 2
        n.1 <- mean(do.call(predict, args)$optimal@Ratio_Stopping_Ind) # saving Ratio Stopping Ind at end of Phase 1

      }
      action[at.risk == 0] <- NA
    } else {
      print("policy is null: generating treatment from rbinom with propensity")
      propensity = predPropensityFn(covariate = covariate)
      action = suppressWarnings(rbinom(N, 1, propensity) * 2 - 1) # 1 for aggressive and 0 for gentle
      # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
      action[at.risk == 0] <- NA
    }

    # hazard
    pred.hazard1 = predHazardFn(action = action, covariate = covariate, cause = 1)
    pred.hazard2 = predHazardFn(action = action, covariate = covariate, cause = 2)

    # generating failure time from cause 1 based on Fine-Gray paper (1999):
    # print(N)
    u1 <- runif(N)  # Observed value
    # print(u1)
    failure_t1 = failure_t2 = obs_time_failureCR = numeric(0)
    for (i in 1:N){
      u1i = u1[i]
      pred.hazard1i = pred.hazard1[i]
      # Use backsolve_t function
      failure_t1[i] <- backsolve_t1(u1i, mass_p, pred.hazard1i)
    }
    failure_t2_rate = exp(pred.hazard2)
    failure_t2 = rexp(N, rate = failure_t2_rate)


    ## summary statistics
    X.cause = pmin(failure_t1, failure_t2)
    censor.time = censor_time
    message(censor.time)
    X = pmin(X.cause, censor.time)
    # print(X.cause)
    D.1 = (failure_t1 <= censor.time & X.cause == failure_t1)
    D.2 = (failure_t2 <= censor.time & X.cause == failure_t2)
    D.0 = 1*D.1+2*D.2
    delta = (X.cause <= censor.time) #overall survival

    # if (full_data == 1) {
    #   return(list(statistics = output,
    #               times = c(failure.time1 = failure_t1, failure.time2 = failure_t2,
    #                         censor.time = censor.time),
    #               D.1 = D.1, D.2 = D.2, D.0 = D.0))
    # } else if (full_data == 2) {
    #   return(c(output, failure.time1 = failure_t1, failure.time2 = failure_t2,
    #            censor.time = censor.time))
    # } else {
    #   return(output)
    # }

    # # tmp3 <<- phase.output
    # print(104)
    # ## bookkeeping and reassigning
    # print(N)
    # print(head(X))
    # print(head(D.0))
    # print(head(action))
    # print(head(at.risk))
    # print(head(output))
    output[, c("event.time", "status",
               # "failure.time1", "failure.time2",
               "action", "at.risk")] <- cbind(X*365, D.0, action, at.risk)
    at.risk[is.na(at.risk)] = 0

    print(output)

