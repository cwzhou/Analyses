#' @param ncov number of covariates
#' @return event.time observed time
#' @return delta 1(min(death, treatment) <= censor)
#' @at.risk Availability at the beginning of the phase
#' @covariates(Z1-Zp) covariate values at the beginning of the phase
#'
Dynamics <-
  function(N = 100,
           u1 = u1,  # Observed value
           u2 = u2,
           u3 = u3,
           n.phases = 2,
           tau = tau,          # structural parameters
           covariate = NULL,
           ncov = 5,
           Sig = diag(ncov) + 0.2 - diag(ncov) * 0.2,      # covariate structure
           mass_p = 0.3,
           censor_time=NULL,
           at.risk = rep(1, N),
           predHazardFn,      # list of predictor functions
           predPropensityFn,  # optimal rule (if !is.null, propensity scores are ignored.) for value calculation
           policy = NULL, # policy is an object of rsf.obj list.
           priority_cause = 1,
           full_data = 0,
           generate_failure_method = NULL,
           crit.eval.os = "mean",
           t0.os = NULL,
           crit.eval.cif = "mean",
           t0.cif = NULL,
           evaluate = FALSE
           # seed1 = seed1
  ) {
    # message("Starting F01.DynamicsCR.R")
    # N = arg.obs$N
    # at.risk = arg.obs$at.risk; ncov = arg.obs$ncov; corr = arg.obs$corr;
    # at.risk = rep(1, N)
    # N = arg.obs$N; n.phases = arg.obs$n.phases; tau = arg.obs$tau;
    # tick = arg.obs$tick; hidden_data = arg.obs$hidden_data; printFlag = arg.obs$printFlag;
    # predHazardFn = arg.obs$predHazardFn; predPropensityFn = arg.obs$predPropensityFn;
    # predCensorFn = arg.obs$predCensorFn; surv.previous = rep(0, N)
    # rho = NULL; omega = NULL; covariate = NULL;
    # Sig = diag(ncov) + 0.2 - diag(ncov) * 0.2
    # censor_time = trunc_cens_time

    # message("N:",N)
    # message("n.phases:",n.phases)
    # message("tau:", tau)
    # message("covariate:", covariate)
    # message("ncov:", ncov)
    # message("mass_p:", mass_p)
    # message("full_data:", full_data)
    # message("censor_time:", censor_time)
    # message("policy:", policy)
    # message("seed1: ", seed1)
    # set.seed(seed1)
    # print(N)
    nphases = 2
    if (length(at.risk) == 1) at.risk = rep(at.risk, N)
    if (length(at.risk) != N)  stop("length of at.risk should match.")
    if (!is.null(covariate)) if (dim(covariate)[1] != N) stop("length of at.risk should match.")

    # initial state
    if (is.null(tau)){
      message("tau is MISSING so we set it to tau = 10")
      tau = 10
    }
    if (is.null(covariate))     covariate = mvrnorm(N, rep(0, ncov), Sig)
    if (is.null(colnames(covariate))) colnames(covariate) = paste0("Z", 1:ncov)
    if (is.null(censor_time)){
      print("censor time was not generated earlier, randomly generating from unif(0,tau) - change if needed")
      censor_time = runif(N,min=0,max=tau)
    }
    if (is.null(u1)){
      message("u1 was null so we are generating - make sure its the same for all methods")
      u1 <<- runif(N)  # Observed value
    }
    if (is.null(u2)){
      message("u2 was null so we are generating - make sure its the same for all methods")
      u2 <<- runif(N)  # Observed value
    }

    # skeleton
    if (evaluate == TRUE){
      output <-
        array(NA, dim = c(N, 4 + 1 + ncov),
              dimnames = list(1:N,
                              c("subj.id", "OS_eval", "CIF_eval","action",
                                "at.risk",
                                colnames(covariate))))
    } else{
      output <-
        array(NA, dim = c(N, 3 + 2 + 2 + ncov), #+1 for failure_t3
              dimnames = list(1:N, c("subj.id", "event.time", "status",
                                     "failure_t1", "failure_t2",#"failure_t3",
                                     "action", "at.risk", colnames(covariate))))
    }

    output[, "subj.id"] <- 1:N
    action = rep(NA, N)           # initialize action vector

    output[, "at.risk"]           <- at.risk
    output[, colnames(covariate)] <- covariate

    ## action                (t=0)
    if (!is.null(policy)) {

      if (is.character(policy) == TRUE){
        if (policy == "test_action0"){
          # print("action = -1 ")
          action = rep(-1, N)
        }
        if (policy == "test_action1"){
          # print("action = 1")
          action = rep(1, N)
        }
      } else{
        if (!is.null(attr(policy, "class"))){

          if (attr(policy, "class") %in% c("ITRSurv")) {
            print("policy is from ITRSurv estimation")
            x = output[as.logical(at.risk),] %>% output2observable(evaluate = TRUE)
            x$Trt = x$A
            # print(1)
            tmp_act = matrix(nrow = N, ncol = nphases+1+1)
            colnames(tmp_act) = c("StopatP1", "P1", "P2", "Final")
            # print(colnames(x))
            # below, we want D.1 b/c priority cause?
            for (phase in 1:2){
              # print(colnames(x))
              # message("F01.DynamicsCR.R: line 103")
              # print(sprintf("Phase: %s", phase))
              if (phase == 1){
                x[, c("obs_time", paste0("D.", phase-1), "Trt", "A")] = 0  # dummy values in order for the package not to drop the NA rows.
              } else{
                x[, c("obs_time", paste0("D.", priority_cause), "Trt", "A")] = 0  # dummy values in order for the package not to drop the NA rows.
              }
              # print(colnames(x))
              x = get_all_vars(policy@phaseResults[[phase]]@model, x)
              # View(x)
              # print(colnames(x))
              pol <<- policy
              xx <<- x
              pp <<- phase
              args <- list(policy, newdata = x, Phase = phase, endPoint = "CR")
              # print(2)
              args_tmp <<- args
              # print("start docalling")
              docalling <- do.call(itrSurv::predict, args)
              # print(3)
              if (phase == 1){
                # print(4)
                docalling1 <<- docalling
                tp_surv <<- policy@params@survivalparam@timePoints
                stopatP1 <<- docalling$optimal@Ratio_Stopping_Ind #stop=1 for P1 and stop=0 for P2
                tmp_act[,1] = stopatP1
                tmp_act[,"P1"] <- docalling1$optimal@optimalTx
              } else{
                # print(4.1)
                docalling2 <<- docalling
                tp_end <<- policy@params@endpointparam@timePoints
                tmp_act[,"P2"] <- docalling2$optimal@optimalTx
              }
            } # end of phase-loop
            tmp_act[, 4] <- ifelse(tmp_act[, "StopatP1"] == 1,
                                   tmp_act[, "P1"],
                                   tmp_act[, "P2"])
            tmp_act_itrsurv <<- tmp_act
            action[at.risk != 0] = tmp_act[,4]
            print(6)
            n.1 <- mean(stopatP1==1) # saving Ratio Stopping Ind at end of Phase 1
            # View(as.data.frame(tmp_act))
            # View(as.data.frame(tmp_act_itrsurv))
          }

          else if (attr(policy, "class") %in% c("DTRSurv")) {
            # single stage CR: n.stages = 1
            n.stages = 1
            print("policy is from DTRSurv estimation (Cho)")
            x = output[as.logical(at.risk),] %>% output2observable(evaluate = TRUE)
            x[, c("obs_time", "D.0")] = 0  # dummy values in order for the package not to drop the NA rows.
            if (n.stages == 1){stage = 1}
            x = get_all_vars(policy@stageResults[[stage]]@model, x)
            if (policy_csk@call$pooled == TRUE){
              x = cbind(x, Trt = NA) # i added this myself to get pooled to work
            }
            # View(x)
            args <- list(policy, newdata = x, stage = stage)
            action[at.risk != 0] <- do.call(dtrSurv::predict, args)$optimal@optimalTx
            tmp_act_csk <<- action
          }
          
          action[at.risk == 0] <- NA
        } else{# policy has attr class
          if (is.list(policy)){
            if (policy[[1]] == "pmcr"){
              # PMCR policy
              print("policy is from PMCR estimation")
              # linear decision rule: d(Z) = I(\beta^TZ>0)
              cov2 <<- covariate %>%
                as.data.frame() %>%
                mutate(int = rep(1, nrow(covariate))) %>%
                dplyr::select("int", everything()) %>%
                as.matrix()
              opt.PMCR =
                ifelse((cov2 %*% as.matrix(policy[[2]])) > 0, 1, -1) %>%
                as.vector()
              action = opt.PMCR
            } else if (policy[[1]] == "aipwe"){
              # AIPWE policy
              print("policy is from AIPWE estimation")
              # decision rule: d(Z) = I(\eta_0 + \eta_1^TZ < 0)
              cov2 <<- covariate %>%
                as.matrix()
              eta0 = policy[[2]][1] #
              eta = policy[[2]][-1]
              opt.AIPWE =
                ifelse((eta0 + cov2 %*% as.matrix(eta)) < 0, 1, -1) %>%
                as.vector()
              action = opt.AIPWE

            } else{
              stop("no policy for that method exists. (F01.DynamicsCR.R LINE 209")
            }
          }
        } # policy does not have attr class
        } # policy is not character
      } # null policy
    else {
      print("policy is null: generating treatment from rbinom with propensity")
      propensity = predPropensityFn(covariate = covariate)
      
      if (revision == 1){
        message("CSK in revision setting can't handle when only 1 patient a group")
        repeat{
          action <- suppressWarnings(rbinom(N, 1, propensity) * 2 - 1)
          action[at.risk == 0] <- NA
          
          # check counts, excluding NAs
          if (sum(action == 1, na.rm = TRUE) >= 2 &&
              sum(action == -1, na.rm = TRUE) >= 2) {
            break  # accept and exit loop
          }
        }
        
        print(sum(action == -1, na.rm = TRUE))
        print(sum(action == 1, na.rm = TRUE))
        stop
      } else{
      action = suppressWarnings(rbinom(N, 1, propensity) * 2 - 1) # 1 for aggressive and -1 for gentle
      # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
      action[at.risk == 0] <- NA
      # print(0)
      }
    }

    # we HAVE OUR ACTION!
    action_tmp1 <<- action
    covariate_tmp1 <<- covariate
    u1_tmp1 <<- u1
    u2_tmp1 <<- u2

    # hazards!!!!!!!
    #hazard linear predictor
    pred.hazard1 = predHazardFn(action = action, covariate = covariate, cause = 1)
    # View(cbind(1,covariate))
    # View(action)
    pred.hazard2 = predHazardFn(action = action, covariate = covariate, cause = 2)
    # pred.hazard3 = predHazardFn(action = action, covariate = covariate, cause = 3)
    # View(pred.hazard1)
    # View(pred.hazard2)

    ph1 <<- pred.hazard1
    ph2 <<- pred.hazard2
    ph <<- cbind(ph1, ph2)
    # print(1)
    failure_os = failure_t1 = failure_t2 = failure_t3 = obs_time_failureCR = numeric(0)

    if (generate_failure_method == "fine_gray"){
      message("fine-gray")
      # generating failure time from cause 1 based on Fine-Gray paper (1999):
      # u11<<-u1
      # pred.hazard11 <<- pred.hazard1
      # massp <<- mass_p
      for (i in 1:N){
        # print(i)
        u1i = u1[i]
        pred.hazard1i = pred.hazard1[i]
        # Use backsolve_t function
        failure_t1[i] <- backsolve_t1(u1i, mass_p, pred.hazard1i)
      }
      # print(1)
      # failure_t1 <- failure_t1
      rate1 = exp(pred.hazard1)
      # print(2)
      #failure t2
      constant_t2 = 1 #0.2 #to make failure time 2 smaller (more of it)
      rate2 <- constant_t2 * exp(pred.hazard2)
      # failure_t2 <- rexp(N, rate = failure_t2_rate)*365 # old one before fixing u2
      failure_t2 <- -(1/rate2)*log(1-u2)
      # print(3)
    } else if (generate_failure_method == "simple_exp"){
      message("simple exponential method")
      # Generate 2 exponential RV which depend on covariates,
      # one for the competing event of interest $E_1$
      # and one for other causes, $E_2$.If $E_1<E_2$,
      # then the individual dies of the event of interest at time $E_1$;
      # otherwise the person dies from the other event at time $E_2$.
      # Hopefully get more direct control over the marginal cumulative incidence curves.
      # rate0 <<- exp(pred.hazard1) + exp(pred.hazard2)
      f1_constant = 1#1.2 # higher f1_constant to make failure time 1 lower (more of it)
      rate1 <- f1_constant * exp(pred.hazard1)
      f2_constant = 1#0.2 # lower f2_constant to make failure time 2 larger (less of it)
      rate2 <- f2_constant * exp(pred.hazard2)
      # set.seed(1);failure_os <- rexp(length(rate0), rate = rate0)
      failure_t1 <- -(1/rate1)*log(1-u1)
      failure_t2 <- -(1/rate2)*log(1-u2)
      # print(head(failure_t1))
      } else{
        stop("no method for generating failure times is specified")
      }
    failure_t1 = failure_t1
    failure_t2 = failure_t2
    # print("failure1 and failure2")
    # print(head(cbind(failure_t1, failure_t2)))
    ftdf <<- cbind(failure_t1, failure_t2)
    rates <<- cbind(rate1 = rate1, rate2 = rate2) #rate0 = rate0,
    # print(1)
    ## evaluate ##
    if (evaluate == TRUE){
      # View(rates)
      # View(rate1)
      # View(rate2)
        OS_eval = eval_sims(generate_failure_method = generate_failure_method,
                            eval_ep = "OS", priority_cause = priority_cause,
                            tau = tau,
                            crit.eval = crit.eval.os, t0 = t0.os,
                            rate1 = rate1,
                            rate2 = rate2,
                            mass_p = mass_p)
        # print("OS_eval")
        # print(head(OS_eval))
        CIF_eval = eval_sims(generate_failure_method = generate_failure_method,
                             eval_ep = "CIF", priority_cause = priority_cause,
                             tau = tau,
                             crit.eval = crit.eval.cif, t0 = t0.cif,
                             rate1 = rate1,
                             rate2 = rate2,
                             mass_p = mass_p)
        # print("CIF_eval")
        # print(head(CIF_eval))
    }

    ## summary statistics
    # ft0 <<- failure_os;
    ft1 <<- failure_t1; ft2 <<- failure_t2
    # View(cbind(ft0 = ft0, ft1 = ft1, ft2 = ft2))
    X.cause = pmin(failure_t1, failure_t2)#, failure_t3)
    # Xcause1 <<- failure_t1; Xcause2 <<- failure_t2
    # Xcause3 <<- failure_t3
    censor.time <<- censor_time # this is already truncated by tau
    # print(3)
    X = pmin(X.cause, censor.time)
    D.1 = (failure_t1 <= censor.time & X.cause == failure_t1)
    D.2 = (failure_t2 <= censor.time & X.cause == failure_t2)
    status = 1*D.1+2*D.2
    # print(98)
    # tmp3 <<- phase.output
    # print(99)
    # ## bookkeeping and reassigning
    # print(100)

    if (evaluate == TRUE){
      output[, c("OS_eval","CIF_eval","action")] <-
        cbind(OS_eval, CIF_eval, action)
    } else{
      output[, c("event.time", "status",
                 "failure_t1", "failure_t2",
                 "action")] <- cbind(X, status,
                                     failure_t1, failure_t2,
                                     action)
    }

    # at.risk[is.na(at.risk)] = 0
    # print(101)

    # } # for (1 in 1:nphases)
    attr(output, "tau") = tau
    attr(output, "ncov") = ncov

    if (full_data == 1) {
      # print("full_data == 1")

      f_times <<- cbind(failure.time1 = failure_t1,
                      failure.time2 = failure_t2,
                      censor.time = censor.time,
                      tau = rep(tau, times = length(failure_t1)),
                      status = rep(NA, times = length(failure_t1))) %>% as.data.frame()
      f_times$observed_time = apply(as.matrix(f_times), 1, min)
      f_times$status = status
      output <<- output
      f1 <<- list(statistics = output,
                  times = f_times,
                  status = status,
                  D.1 = D.1,
                  D.2 = D.2)
      return(f1)
    } else if (full_data == 2) {
      return(c(output, failure.time1 = failure_t1, failure.time2 = failure_t2,
               censor.time = censor.time))
    } else {
      return(output)
    }
  }

# backsolving for failure time from cause 1 based off of Fine-Gray 1999 paper
backsolve_t1 <- function(u1i, mass_p, pred.hazard1) {
  # Function to solve for t
  solve_t <- function(t, u1i, mass_p, pred.hazard1) {
    ((1 - ((1 - mass_p * (1 - exp(-t))) ^ exp(pred.hazard1))) / (1 - ((1 - mass_p) ^ exp(pred.hazard1)))) - u1i
  }
  # Specify the interval
  lower_boundary <- 1e-250
  upper_boundary <- 1e250
  # Evaluate the function at the interval's endpoints
  value_at_lower <- solve_t(lower_boundary, u1i, mass_p, pred.hazard1)
  value_at_upper <- solve_t(upper_boundary, u1i, mass_p, pred.hazard1)
  if (is.na(value_at_lower)){
    # Print the values
    # message("u1i:",u1i)
    # message("pred.hazard1:",pred.hazard1)
    }
  # cat("Value at lower boundary:", value_at_lower, "\n")
  # cat("Value at upper boundary:", value_at_upper, "\n")

  # Use uniroot to solve for t
  result1 <- uniroot(solve_t, lower = lower_boundary, upper = upper_boundary, tol = 1e-9, u1i = u1i, mass_p = mass_p, pred.hazard1 = pred.hazard1)
  return(result1$root)
}

# define simulation value formulas for evaluating based on sim setting
# overall survival for evaluation for sims
OS_func <- function(generate_failure_method, t, cause1_prob, lambda1, lambda2) {
  if (generate_failure_method == "simple_exp"){
    # lambda1 = exp(Z %*% beta1.treatmentA)
    # lambda2 = exp(Z %*% beta2.treatmentA)
    # b/c assuming independence between T1 and T2
    rate = lambda1 + lambda2
    St = exp(-rate * t)
  } else if (generate_failure_method == "fine_gray"){
    # print("fg os")
    term1 = (1-cause1_prob*(1-exp(-t)))^lambda1
    term2 = (1-exp(-t*lambda2))*((1-cause1_prob)^lambda1)
    St = term1 - term2
  } else{
    stop("generate_failure_method must be simple_exp or fine_gray.")
  }
  # message("St = ", St)
  return(St)
}
# PC subdistribution for evaluation for sims
subdist_func <- function(generate_failure_method, priority_cause, t, cause1_prob, lambda1, lambda2){
  if (generate_failure_method == "simple_exp"){
    if (priority_cause == 1){
      lambda_pc = lambda1
      lambda_npc = lambda2
    } else{
      lambda_pc = lambda2
      lambda_npc = lambda1
    }
    # Pr(T\leqt, \epsilon=priority_cause): PC subdistribution
    Lambda_sum = lambda_pc + lambda_npc
    term1 = lambda_pc/Lambda_sum
    term2 = -exp(-lambda_npc * t)
    term3 = (lambda_npc/Lambda_sum) * exp(-t * Lambda_sum)
    # View(term1)
    # View(term2)
    # View(term3)
    subdistr = term1 + term2 + term3
  } else if (generate_failure_method == "fine_gray"){
    # print("fg cif")
    if (priority_cause == 1){
      # Pr(T\leq t,\epsilon=1) = 1 - [1- p(1-exp(-t))]^(lambda1)
      subdistr = 1 - (1- cause1_prob*(1-exp(-t)))^(lambda1)
    } else if (priority_cause == 2){
      # Pr(T\leq t,\epsilon=2|Z) = Pr(T\leq t|\epsilon=2,Z)Pr(\epsilon=2|Z)
      # Pr(T\leq t|\epsilon=2,Z)
      cdf_exp = 1-exp(-lambda2*t)
      # Pr(\epsilon = 2 | Z)
      Pr_eps2 = (1-cause1_prob)^lambda1
      subdistr = cdf_exp * Pr_eps2
      #lambda1 = exp(Z%*%beta1.treatmentA)
    }
  } else{
    stop("generate_failure_method must be simple_exp or fine_gray.")
  }
  # message("CIF = ", subdistr)
  return(subdistr)
}

# evaluate
eval_sims = function(generate_failure_method,
                     eval_ep,
                     priority_cause,
                     tau,
                     crit.eval,
                     t0,
                     rate1, rate2,
                     mass_p){
  # message('EVALUATING: DynamicsCR.r line 263')

  if (eval_ep == "OS"){
    # overall survival
    if (crit.eval == "prob"){
      # message("crit.eval is prob")
      if (is.null(t0) | is.na(t0)){
        stop("t0.os needs to be specified for current crit.eval.os in F01 LINE 273")
      }
      OS_eval = OS_func(generate_failure_method, t0, mass_p, rate1, rate2)
    } else{
      # Set the range for integration (e.g., from 0 to infinity)
      OS_result1 = numeric(0)
      for (i in 1:length(rate1)){
        OS_result <- integrate(OS_func, lower = 0, upper = tau,
                               generate_failure_method = generate_failure_method,
                               cause1_prob = mass_p,
                               lambda1 = rate1[i], lambda2 = rate2[i])
        if (i == 1){
          # print("OS testing i = 1")
          # print(OS_result$value)
        }
        OS_result1[i] = OS_result$value
        if (i == 1){
        #   message("rate1:",rate1[i])
        #   message("rate2:",rate2[i])
          # print(OS_func)
          # message("OS_result$value:",OS_result$value)
        }
      }
      OS_eval = OS_result1
      testt <<- OS_func(generate_failure_method, t0, mass_p, rate1, rate2)
      # View(testt)
    }
    eval_result = OS_eval
  } else if (eval_ep == "CIF"){
    # endpoint
    if (crit.eval == "prob"){
      # message("crit.eval is prob")
      if (is.null(t0) | is.na(t0)){
        stop("t0.cif needs to be specified for current crit.eval.cif in F01 LINE 273")
      }
      CIF_eval = subdist_func(generate_failure_method, priority_cause, t0, mass_p, rate1, rate2)
    } else{
      # Set the range for integration (e.g., from 0 to infinity)
      CIF_result1 = numeric(0)
      for (i in 1:length(rate1)){
        CIF_result <- integrate(subdist_func, lower = 0, upper = tau,
                                generate_failure_method = generate_failure_method,
                                priority_cause = priority_cause,
                                cause1_prob = mass_p,
                                lambda1 = rate1[i], lambda2 = rate2[i])
        if (i == 1){
          # print("CIF testing i = 1")
          # print(CIF_result$value)
          # stop()
        }
        CIF_result1[i] = CIF_result$value

      }
      CIF_eval = CIF_result1
      testtcr <<- subdist_func(generate_failure_method, priority_cause, t0, mass_p, rate1, rate2)
      # View(testtcr)
    }
    eval_result = CIF_eval
    } else{
    stop("eval_ep needs to be one of OS or CIF")
    }
  return(eval_result)
}


# message("End of F01.DynamicsCR.R")


# End of script -------------------------------------------------------------

