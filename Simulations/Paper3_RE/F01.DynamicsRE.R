#' @param ncov number of covariates
#' @return event.time.start observed time
#' @return event.time.end observed time
#' @return indR: 1(re <= min(death,censor))
#' @return indD: 1(death <= censor)
#' We assume death and RE cannot occur at same time.
#' @covariates(Z1-Zp) covariate values at the beginning of the phase
#'
Dynamics <-
  function(N = 100, #number of SUBJECTS
           u1 = u1,
           u2 = u2,
           u3 = u3,
           tau = tau0,          # structural parameters
           covariate = NULL,
           ncov = 5,
           Sig = diag(ncov) + 0.2 - diag(ncov) * 0.2,      # covariate structure
           # gapparam1, # alpha1/2 depends on treatment, so moved gapparam1/2 to within F01.DynamicsRE.R, as of Jan 3, 2025
           # gapparam2,
           G = G,
           censor_time=NULL,
           predHazardFn_D,      # list of predictor functions
           predHazardFn_R,
           predPropensityFn,  # optimal rule (if !is.null, propensity scores are ignored.) for value calculation
           policy = NULL, # policy is an object of rsf.obj list.
           crit.eval.surv = "mean",
           t0.surv = NULL,
           crit.eval.endpoint = "mean",
           t0.endpoint = NULL,
           evaluate = FALSE
  ) {
    # message("Starting F01.DynamicsRE.R")
    # message("N:",N)
    # message("tau:", tau)
    # message("covariate:", covariate)
    # message("ncov:", ncov)
    # message("full_data:", full_data)
    # message("censor_time:", censor_time)
    # message("policy:", policy)

    if (!is.null(covariate)){
        if (dim(covariate)[1] != N){
          print(dim(covariate)[1])
          print(N)
          stop("length of cov should match.")
        }
    }

    # initial state
    if (is.null(tau)){
      message("tau is NULL so we arbitraily set it to tau = 10")
      tau = 10
    }
    if (is.null(covariate)){
      message("covariate is NULL so we arbitraily set it to covariate = mvrnorm(N, rep(0, ncov), Sig)")
      covariate = mvrnorm(N, rep(0, ncov), Sig)
    }
    if (is.null(colnames(covariate)) & ncov > 1) colnames(covariate) = paste0("Z", 1:ncov)
    if (is.null(censor_time)){
      print("censor_time is NULL, randomly generating from unif(0,tau) - change if needed")
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

    epName1 = "IndR"
    # skeleton
    if (evaluate == TRUE){
      output <-
        array(NA, dim = c(N, 7 + ncov),
              dimnames = list(1:N, c("subj.id",
                                     "L_open", "R_closed",
                                     "IndD", epName1,
                                     "obs.failure",
                                     "action", colnames(covariate))))
      # output <-
      #   array(NA, dim = c(N, 4 + ncov),
      #         dimnames = list(1:N,
      #                         c("subj.id", "Surv_eval", "Endpoint_eval","action",
      #                           colnames(covariate))))
    } else{
      output <-
        array(NA, dim = c(N, 7 + ncov),
              dimnames = list(1:N, c("subj.id",
                                     "L_open", "R_closed",
                                     "IndD", epName1,
                                     "obs.failure",
                                     "action", colnames(covariate))))
    }

    output[, "subj.id"] <- 1:N
    action = rep(NA, N)           # initialize action vector
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
            op <<- output
            x = output %>% as.data.frame() #%>% output2observable(evaluate = TRUE)
            x$Trt = x$A
            # print(1)
            tmp_act <<- matrix(nrow = N, ncol = 2+1+1)
            colnames(tmp_act) = c("StopatP1", "P1", "P2", "Final")
            # below, we want D.1 b/c priority cause?
            for (phase in 1:2){
              message("F01.DynamicsRE.R: line 116 for Phase ", phase)
              x[, c("L_open","R_closed",
                    "IndD", epName1, "obs.failure",
                    "action")] = 0  # dummy values in order for the package not to drop the NA rows.
              # print(colnames(x))
              mm <<- policy@phaseResults[[phase]]@model
              xx <<- x
              x = get_all_vars(policy@phaseResults[[phase]]@model, x)
              # print("col2")
              # print(colnames(x))
              poll <<- policy
              phasee <<- phase
              neww <<- x
              # print(2)
              args <- c(policy,list(newdata = x, Phase = phase, epName1 = epName1, endPoint = endpoint))
              args_tmp <<- args
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
            action = tmp_act[,4]
            # print(6)
            n.1 <- mean(stopatP1==1) # saving Ratio Stopping Ind at end of Phase 1
            # View(as.data.frame(tmp_act))

        } else{# policy has attr class
            stop("no policy for that method exists. (F01.DynamicsRE.R LINE 158")
          }
        } # policy does not have attr class
      } # policy is not character
    } # null policy
    else {
      print("policy is null: generating treatment from rbinom with propensity")
      # print(predPropensityFn)
      propensity = predPropensityFn(covariate = covariate)
      # print(propensity)
      action = suppressWarnings(rbinom(N, 1, propensity)) # 1 for aggressive and 0 for gentle
      # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
      # print(0)
      prop_tmp1 <<- propensity
    }

    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################
    # we HAVE OUR ACTION!
    action_tmp1 <<- action
    covariate_tmp1 <<- covariate
    # u1_tmp1 <<- u1
    # u2_tmp1 <<- u2

    # generating gap time 1 using Gumbel bivariate exponential model

    # hazards!
    # lambda_D = lambda_0D*exp(t(beta_D)%*%z + omega_D * A + A * t(gamma_D) %*% z) #no intercept for z b/c survival
    pred.hazard1 <<- predHazardFn_D(action = action, covariate = covariate)
    # print(pred.hazard1)
    # lambda_R = lambda_0R*exp(t(beta_R)%*%z + omega_R * A + A * t(gamma_R) %*% z)
    pred.hazard2 <<- predHazardFn_R(action = action, covariate = covariate)
    # print(pred.hazard2)
    # message('action', action)
    # message('covariate', covariate)

    if (evaluate == TRUE){
      # message("original G was: ", G)
      G = G #999
      # message("But now setting G to be really high: G = ", G)
    } else{
      G = G
    }

    #if alpha=0 then independent and don't need to do this.
    gapparam1 <<- ifelse(action == 1, 0.8, 0.5) #0.5, 0.1) #rho for failure
    # print(gapparam1)
    gapparam2 <<- ifelse(action == 1, 0.2, 0.7) #0.4, 0.7) #rho for gap times
    alpha1 = 4*gapparam1 # alpha1/2 depends on treatment, so moved gapparam1/2 to F01.DynamicsRE.R, as of Jan 3, 2025
    alpha2 = 4*gapparam2

    if (any(alpha1 == 0)) {
      print("At least one subject has alpha1 = 0 so creating rateD")
      rateD = predHazardFn_D(action = 0, covariate = covariate)
    } else{
      rateD = NULL
    }
    if (any(alpha2 == 0)) {
      print("At least one subject has alpha2 = 0 so creating rateR")
      rateR = predHazardFn_R(action = 0, covariate = covariate)
    } else{
      rateR = NULL
    }
    # print("============== Gumbel Bivariate Exponential ==============")
    # generate first gap time using marginal exp(lamR) distribution
    # gaptime1 = rexp(N,lambda_R)
    # gaptime1 <<- rexp(N,pred.hazard2)
    gaptime1 <<- -(1/pred.hazard2)*log(1-u1)
    # print(sprintf("gaptime1: %s", head(gaptime1)))

    # print("Generating Failure Time")
    # generate conditional failure time (given gaptime1)
    tt_fail <<- as.numeric(GumbelBiExp(N=N,lambda_D=pred.hazard1,lambda_R=pred.hazard2,
                                       rate = rateD, # needed for if alpha1[i] is 0
                                       alpha=alpha1,y_type=1,y=gaptime1,u=u2)$tt_fail)
    # print(sprintf("failure time: %s", head(tt_fail)))

    # Generate G conditional gaptime values for each subject
    if (G>1){
      # gaptimes <- generate_gaptime(N, lambda_D, lambda_R, alpha2, G)
      gaptimes <- generate_gaptime(N, pred.hazard1, pred.hazard2, rate = rateR, alpha2, G, u = u3)
      gaptimes[,1] = gaptime1
      gap_names <- paste0("gaptime", 1:ncol(gaptimes))
      colnames(gaptimes) <- gap_names
    } else if (G==1){
      gaptimes = gaptime1 %>% as.data.frame()
      # print(gaptimes)
      colnames(gaptimes) = "gaptime1"
    }

    # Create a dataframe with the gaptime values
    gap_df <- data.frame(ID = 1:N, gaptime0 = rep(0,N),gaptimes);#print(head(gap_df))
    # Pivot the dataset to a longitudinal format
    gap_df_longitudinal <- gap_df %>%
      pivot_longer(cols = starts_with("gaptime"), names_to = "Gap", values_to = "Time_Gap") %>%
      mutate(Gap = as.numeric(gsub("gaptime", "", Gap)));
    # print(head(gap_df_longitudinal))

    # Calculate the recurrent times using cumulative sum per ID
    recurrent_df_longitudinal <- gap_df_longitudinal %>%
      group_by(ID) %>%
      mutate(Time_Recurrent = cumsum(Time_Gap),
             Recurrent = Gap) %>%
      ungroup() %>%
      filter(Gap != 0) %>%
      as.data.frame() %>%
      mutate(Label = paste0("Recurrent", Recurrent)) %>%
      dplyr::select(-c(Gap, Time_Gap, Recurrent)) %>%
      rename(Time = Time_Recurrent)

    # # generate conditional gaptime2 (given gaptime1)
    # gaptime2 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime1,u=u3)$tt)
    # gaptime3 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime2,u=u4)$tt)
    # gaptime4 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime3,u=u5)$tt)
    # gap_df = data.frame(ID = c(1:N), gaptime1 = gaptime1, gaptime2 = gaptime2, gaptime3 = gaptime3, gaptime4 = gaptime4)

    # # print(alpha1)
    # if (alpha1 == 0){
    #   print("Independence")
    #   # tt_fail = rexp(N,lambda_0D*exp(t(beta_D)*z))
    #   # tt_fail <<- rexp(N, predHazardFn_D(action = 0, covariate = covariate))
    #   rateD = predHazardFn_D(action = 0, covariate = covariate)
    #   tt_fail <<- -(1/rateD)*log(1-u1)
    # } else{
    #   # print("============== Gumbel Bivariate Exponential ==============")
    #
    #   # generate first gap time using marginal exp(lamR) distribution
    #   # gaptime1 = rexp(N,lambda_R)
    #   # gaptime1 <<- rexp(N,pred.hazard2)
    #   gaptime1 <<- -(1/pred.hazard2)*log(1-u1)
    #   # print(sprintf("gaptime1: %s", head(gaptime1)))
    #
    #   # print("Generating Failure Time")
    #   # generate conditional failure time (given gaptime1)
    #   tt_fail <<- as.numeric(GumbelBiExp(N=N,lambda_D=pred.hazard1,lambda_R=pred.hazard2,
    #                                      alpha=alpha1,y_type=1,y=gaptime1,u=u2)$tt)
    #   # print(sprintf("failure time: %s", head(tt_fail)))
    #
    #   # print("Checking Plot")
    #   # plot_check = Check_GumbelBiExp(N=N,tt=tt_fail);
    #   # print(plot_check)
    #
    #   # Generate G conditional gaptime values for each subject
    #   if (G>1){
    #     # gaptimes <- generate_gaptime(N, lambda_D, lambda_R, alpha2, G)
    #     gaptimes <- generate_gaptime(N, pred.hazard1, pred.hazard2, alpha2, G, u = u3)
    #     gaptimes[,1] = gaptime1
    #     gap_names <- paste0("gaptime", 1:ncol(gaptimes))
    #     colnames(gaptimes) <- gap_names
    #   } else if (G==1){
    #     gaptimes = gaptime1 %>% as.data.frame()
    #     # print(gaptimes)
    #     colnames(gaptimes) = "gaptime1"
    #   }
    #
    #   # Create a dataframe with the gaptime values
    #   gap_df <- data.frame(ID = 1:N, gaptime0 = rep(0,N),gaptimes);#print(head(gap_df))
    #   # Pivot the dataset to a longitudinal format
    #   gap_df_longitudinal <- gap_df %>%
    #     pivot_longer(cols = starts_with("gaptime"), names_to = "Gap", values_to = "Time_Gap") %>%
    #     mutate(Gap = as.numeric(gsub("gaptime", "", Gap)));
    #   # print(head(gap_df_longitudinal))
    #
    #   # Calculate the recurrent times using cumulative sum per ID
    #   recurrent_df_longitudinal <- gap_df_longitudinal %>%
    #     group_by(ID) %>%
    #     mutate(Time_Recurrent = cumsum(Time_Gap),
    #            Recurrent = Gap) %>%
    #     ungroup() %>%
    #     filter(Gap != 0) %>%
    #     as.data.frame() %>%
    #     mutate(Label = paste0("Recurrent", Recurrent)) %>%
    #     dplyr::select(-c(Gap, Time_Gap, Recurrent)) %>%
    #     rename(Time = Time_Recurrent)
    #
    #   # # generate conditional gaptime2 (given gaptime1)
    #   # gaptime2 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime1,u=u3)$tt)
    #   # gaptime3 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime2,u=u4)$tt)
    #   # gaptime4 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime3,u=u5)$tt)
    #   # gap_df = data.frame(ID = c(1:N), gaptime1 = gaptime1, gaptime2 = gaptime2, gaptime3 = gaptime3, gaptime4 = gaptime4)
    # } # end of Bivariate Gumbel

    # if (evaluate == TRUE){
    #   Surv_eval = eval_sims(eval_ep = "Surv",
    #                       tau = tau,
    #                       crit.eval = crit.eval.surv,
    #                       t0 = t0.surv,
    #                       rate1 = rate1,
    #                       rate2 = rate2)
    #   print("Surv_eval")
    #   print(head(Surv_eval))
    #   Endpoint_eval = eval_sims(eval_ep = "MFF",
    #                        tau = tau,
    #                        crit.eval = crit.eval.endpoint,
    #                        t0 = t0.endpoint,
    #                        rate1 = rate1,
    #                        rate2 = rate2)
    #   print("Endpoint_eval")
    #   print(head(Endpoint_eval))
    # }

    #steps to put together in one dataset
    trunc_cens_time = pmin(censor_time, tau0)
    df_times <<- data.frame(ID = c(1:N), Time_Failure=tt_fail, Time_Censor=censor_time, Time_Tau=tau0, Truncated_Censor_Time = trunc_cens_time)
    # View(df_times)

    # survival dataset (1row/person)
    df_cov = data.frame(ID = c(1:N), covariate, Trt = action)
    # View(df_cov)
    df_surv1 = df_times %>% mutate(obs_time = pmin(Time_Failure, Time_Censor, Time_Tau),
                                   indD = ifelse(Time_Failure <= pmin(Time_Censor, Time_Tau), 1, 0))
    times_act <<- inner_join(df_cov, df_surv1, by = "ID")

    # View(df_surv1)
    df_surv2 = inner_join(df_cov, df_surv1, by = "ID") %>%
      dplyr::select(-c(Time_Failure, Time_Censor, Time_Tau))
      # dplyr::select(ID, all_of(covariate), Trt, obs_time, indD)
    # View(df_surv2)

    # recurrent dataset
    df_times_long <<- df_times %>%
      dplyr::select(-c(Truncated_Censor_Time)) %>%
      pivot_longer(cols = starts_with("Time_"), names_to = "Label", values_to = "Time") %>%
      mutate(Label = gsub("Time_", "", Label)) %>%
      as.data.frame();
    # View(df_times_long)
    df_long <<- rbind(recurrent_df_longitudinal, df_times_long) %>%
      dplyr::select(ID, Label, Time) %>%
      group_by(ID) %>%
      arrange(ID, Time)
    # View(df_long)

    # find min(failure, censoring, tau0)
    df_min <<- df_long %>%
      group_by(ID) %>%
      mutate(failure_status = ifelse(Label == "Failure", 1, ifelse(Label == "Censor" | Label == "Tau", 0, 99))) %>%
      filter(Label %in% c("Censor", "Failure", "Tau")) %>%
      summarise(min = min(Time),
                failure_status_raw = failure_status[which.min(Time)],
                Label = Label[which.min(Time)])
    # View(df_min)

    # combine min with long dataset
    dataset0 = inner_join(df_long, df_min, by = "ID") %>%
      # only keep observed data
      filter(Time <= min) %>%
      as.data.frame() %>%
      dplyr::select(-Label.y) %>%
      group_by(ID) %>%
      # create dummy var for first row within ID, recurrence, failure, censoring
      mutate(FirstRow = ifelse(row_number() == 1, 1, 0),
             contains_recurrent = (grepl("recurrent", Label.x, ignore.case = TRUE)),
             contains_failure = (grepl("failure", Label.x, ignore.case = TRUE)),
             contains_censor = (grepl("censor", Label.x, ignore.case = TRUE))) %>%
      # use lag to move future times to "L_open"
      mutate(L_open = ifelse(FirstRow ==1, 0, lag(Time)), #open parentheses ("(")
             R_closed = Time, #closed brack ("]")
             IndR = ifelse(contains_recurrent == TRUE, 1, 0),
             IndD = ifelse(contains_failure == TRUE, 1, 0)) %>%
      rename(R_Label = Label.x) %>%
      dplyr::select(ID, L_open, R_closed, IndR, IndD, R_Label) %>%
      as.data.frame()
    # View(dataset0)

    #add in covariates
    dataset = merge(dataset0, df_cov, by = "ID") #%>%
      # dplyr::select(ID, L_open, R_closed, IndR, IndD, covariate, Trt, R_Label)
    # View(dataset)

    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    if (evaluate == TRUE){
      # output[, c("Surv_eval","Endpoint_eval","action")] <-
      #   cbind(Surv_eval, Endpoint_eval, action)
      # attr(output, "tau") = tau
      # attr(output, "ncov") = ncov
      # return(output)

      today_date <- format(Sys.Date(), "%Y-%m-%d")
      name = sprintf("Evaulation_Dataset_%s_N%s_G%s_tau%s_ncov%s",
                     today_date,
                     n, G, tau0, ncov)
      return(list(dataset_recurrent=dataset,dataset_survival=df_surv2,name=name))


    } else{
      today_date <- format(Sys.Date(), "%Y-%m-%d")
      name = sprintf("Dataset_%s_N%s_G%s_tau%s",
                     today_date,
                     n, G, tau0)
      return(list(dataset_recurrent=dataset,dataset_survival=df_surv2,name=name))
      # output[, c("L_open", "R_closed",
      #            "IndD", epName1,
      #            "obs.failure",
      #            "action")] <- cbind(X.start, X.end,
      #                                IndD, IndRE,
      #                                failure_t1, failure_t2,
      #                                action)
    }

  } # end of Dynamics function

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# below are functions used in Dynamics() function
# evaluate
eval_sims = function(eval_ep,
                     tau,
                     crit.eval,
                     t0,
                     rate1,
                     rate2){
  message('EVALUATING: F01.DynamicsRE.R, line 389')

  if (eval_ep == "Surv"){
    # survival
    if (crit.eval == "prob"){
      # message("crit.eval is prob")
      if (is.null(t0) | is.na(t0)){
        stop("t0.surv needs to be specified for current crit.eval.surv in F01 LINE 273")
      }
      Surv_eval = Surv_func(t0, rate1)
    } else{
      # Set the range for integration
      Surv_result1 = numeric(0)
      for (i in 1:length(rate1)){
        Surv_result <- integrate(Surv_func, lower = 0, upper = tau,
                               lambda1 = rate1[i])
        Surv_result1[i] = Surv_result$value
      }
      Surv_eval = Surv_result1
      testt <<- Surv_func(t0, rate1)
    }
    eval_result = Surv_eval
  } else if (eval_ep == "MFF"){
    # endpoint
    if (crit.eval == "prob"){
      # message("crit.eval is prob")
      if (is.null(t0) | is.na(t0)){
        stop("t0.endpoint needs to be specified for current crit.eval.endpoint in F01 LINE 273")
      }
      Endpoint_eval = mff_func(t0, rate1, rate2)
    } else{
      # Set the range for integration
      Endpoint_result1 = numeric(0)
      for (i in 1:length(rate1)){
        # pretty sure we dont want integration here, so fix this. this is for CR which is not right for RE.
        Endpoint_result <- integrate(mff_func, lower = 0, upper = tau,
                                lambda1 = rate1[i], lambda2 = rate2[i])
        Endpoint_result1[i] = Endpoint_result$value
      }
      Endpoint_eval = Endpoint_result1
      testtre <<- mff_func(t0, rate1, rate2)
    }
    eval_result = Endpoint_eval
  } else{
    stop("eval_ep needs to be one of Surv or MFF")
  }
  return(eval_result)
}

# define simulation value formulas for evaluating based on sim setting
# need to fix below, bc rn its for CR
# survival for evaluation for sims
Surv_func <- function(t, lambda1) {
  if (alpha1 == 0){ # independent
    St = exp(-lambda1 * t)
  } else{ # gumbel bivariate conditional survival time
    St = exp(-t)*(1+alpha1-2*alpha1*exp(-lambda1*y)) + alpha1*exp(-2*t)*(1-2*exp(-lambda1*y))
  }

  message("St = ", St)
  return(St)
}

# need to fix below, because rn its subdistr.
# mean freq func for evaluation for sims
mff_func <- function(t, lambda1, lambda2){
  lambda_pc = lambda2
  lambda_npc = lambda1
  # Pr(T\leqt, \epsilon=priority_cause): PC subdistribution
  Lambda_sum = lambda_pc + lambda_npc
  term1 = lambda_pc/Lambda_sum
  term2 = -exp(-lambda_npc * t)
  term3 = (lambda_npc/Lambda_sum) * exp(-t * Lambda_sum)
  mff = term1 + term2 + term3
  message("MFF = ", mff)
  return(mff)
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

message("End of F01.DynamicsRE.R")


# End of script -------------------------------------------------------------
