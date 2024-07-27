# Title: Libraries and Functions for Simulation Scripts
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
} else{
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper1_CR")
}
source("F00.Simulation_Libraries.R")
source("F02.ComparatorMethod_Functions.R")
source("F01.multiPhaseDynamicsCR.R")
# Functions -----------------------------------------------------------------
gdata_CR <- function(N=10,
                     u1 = NULL,
                     u2 = NULL,
                     u3 = NULL,
                     ncov = 5,
                     mass_p = 0.3,
                     predHazardFn, predPropensityFn,
                     ztype=0, zparam=0.5,
                     ctype=0,cparam=2,censor_min=0,censor_max=NULL,
                     num_A=2,tau=5,
                     policy=NULL, priority_cause = 1,
                     full_data=0,
                     generate_failure_method = NULL,
                     crit.eval.os = "mean",
                     t0.os = NULL,
                     crit.eval.cif = "mean",
                     t0.cif = NULL,
                     evaluate = FALSE){
  # N = arg.csk$N; ncov = arg.csk$ncov; mass_p = arg.csk$mass_p;
  # ctype = arg.csk$ctype
  # cparam=2;censor_min=0;censor_max=NULL;
  # ztype = 0; zparam = 0.5; num_A = 2
  # tau = arg.csk$tau; policy = arg.csk$policy; seed1 = arg.csk$seed1
  # message("N:",N)
  # message("ncov:", ncov)
  # message("mass_p:", mass_p)
  # message("ztype:", ztype)
  # message("zparam:", zparam)
  # message("ctype:", ctype)
  # message("cparam:", cparam)
  # message("censor_min:", censor_min)
  # message("censor_max:", censor_max)
  # message("num_A:", num_A)
  # message("tau:", tau)
  # message("seed1: ", seed1)
  if (is.null(u1)){
    print("missing u1 --- make sure this is right...")
    u1 <<- runif(N)  # Observed value
  }
  if (is.null(u2)){
    print("missing u2 --- make sure this is right...")
    u2 <<- runif(N)  # Observed value
  }
  if (is.null(u3)){
    print("missing u3 --- make sure this is right...")
    u3 <<- runif(N)  # Observed value
  }
  # print(beta1.hazard1)
  # ztype indicates the distribution for covariate z: 0=normal(0,1),1=binary(zparam),2=uniform(0,1)
  # ctype indicates the distribution for censoring (0=exponential,1=uniform,
  #                                                 9=no censoring)

  # name = sprintf("CRDataset_N%s_G%s_A%s_lambda0D%sBetaD%somegaD%sgammaR%s_lambda0R%sBetaR%somegaR%sgammaR%s_rho1%s_rho2%s_tau%s",
  #                N, G, num_A,
  #                round(lambda_0D,1), round(beta_D,1), round(omega_D,1), round(gamma_D,1), round(lambda_0R,1), round(beta_R,1), round(omega_R,1), round(gamma_R,1),
  #                gapparam1, gapparam2, tau)
  # df_cov = df_times = df_times0 = df_surv1 = df_surv2 = data.frame(NA)
  # message("dim of df_times0 is first:", dim(df_times0))
  # print("----------------")
  # generating covariates (nxp matrix)
  # print(N)
  if (ztype == 0){
    # print(N)
    # print(ncov)
    print("ztype=normal")
    z <- matrix(rnorm(N * ncov), nrow = N, ncol = ncov) # positive for when testing but regular when we want zom to be not as good
    # standardized_z = scale(z)
    # print(head(z))
    # print(head(standardized_z))
    # z = standardized_z
    # print("standardized z")

  }
  # if (ztype == 0) {
  #   z <- matrix(0, nrow = N, ncol = ncov)
  #   # First covariate is abs(rnorm) (positive)
  #   z[, 1] <- abs(rnorm(N))
  #   # Other covariates are regular rnorm
  #   z[, -1] <- matrix(rnorm((N * (ncov - 1))), nrow = N, ncol = ncov - 1)
  # }
  if (ztype == 1){z <- matrix(rbinom(N*ncov, 1, zparam),nrow=N,ncol=ncov)}
  if (ztype == 2){z <- matrix(runif(N*ncov),nrow=N,ncol=ncov)}
  if (is.null(colnames(z))) colnames(z) = paste0("Z", 1:ncov)

  # generating censoring time
  if (ctype == 0){
    message("censoring: exp")
    cc <- rexp(N,cparam)
    }
  if (ctype == 1){
    message("censoring: unif")
    if(is.null(censor_max)){censor_max = tau}
    if(is.null(censor_min)){censor_min = 0}
    message("censor_max is", censor_max)
    cc <- runif(N,min=censor_min,max=censor_max)
    # View(cc)
    }
  if (ctype == 99){
    print("!!!!! no censoring. !!!!!")
    cc = rep(10^250, N) #arbitrary large number so never censored via cc = censor time (could still be censored by tau though so fix in future)
  }
  trunc_cens_time = pmin(cc, tau)
  # print(trunc_cens_time)
  # print(sprintf("censoring time cc: %s",cc))
  # print(cc)
  # generating treatment
  # A = rbinom(N, num_A-1, 0.5)
  # print(sprintf("treatments A: %s",A))
  # print(N)
  df_multi <<- multiPhaseDynamics(N=N, u1 = u1, u2=u2, u3 = u3,
                         tau=tau, covariate = z, ncov = ncov,
                         mass_p = mass_p,
                         censor_time = trunc_cens_time,
                         predHazardFn = predHazardFn,
                         predPropensityFn = predPropensityFn,
                         policy=policy, priority_cause=priority_cause,
                         full_data = full_data,
                         generate_failure_method = generate_failure_method,
                         crit.eval.os = crit.eval.os,
                         t0.os = t0.os,
                         crit.eval.cif = crit.eval.cif,
                         t0.cif = t0.cif,
                         evaluate = evaluate
                           ) %>% as.data.frame()



  # propensity = predPropensityFn(covariate = z)
  # A = suppressWarnings(rbinom(N, 1, propensity) * 2 - 1) # 1 for aggressive and 0 for gentle
  # pred.hazard1 = predHazardFn(action = A, covariate = z, cause = 1)
  # pred.hazard2 = predHazardFn(action = A, covariate = z, cause = 2)
  #
  # df_cov = data.frame(c(1:N), z, A)
  # colnames(df_cov) = c("ID", colnames(z), "Trt")
  #
  # # generating failure time from cause 1 based on Fine-Gray paper (1999):
  # u1 <- runif(N)  # Observed value
  # failure_t1 = failure_t2 = obs_time_failureCR = numeric(0)
  # for (i in 1:N){
  #   u1i = u1[i]
  #   pred.hazard1i = pred.hazard1[i]
  #   # Use backsolve_t function
  #   failure_t1[i] <- backsolve_t1(u1i, mass_p, pred.hazard1i)
  # }
  # failure_t2_rate = exp(pred.hazard2)
  # failure_t2 = rexp(N, rate = failure_t2_rate)
  # pred.hazard2 <<- pred.hazard2
  # failre_t2_rate <<- failure_t2_rate
  # failre_t2 <<- failure_t2
  # # plot(failure_t2, main = "test: failure from cause 2")
  # # message("N:",N)
  # # print(length(failure_t1));print(length(failure_t2));
  # # print(length(cc));print(length(tau))
  # df_times0 = data.frame(ID = 1:N,
  #                   Time_Failure1 = failure_t1,
  #                   Time_Failure2 = failure_t2,
  #                   Time_Censor = cc,
  #                   Time_Tau = tau
  #                   )
  # df_times0 <<-df_times0
  # # print(length(df_times0$Time_Failure2));print(length(df_times0$Time_Failure1));
  # # print(length(df_times0$Time_Censor));print(length(df_times0$Time_Tau))
  # df_times = df_times0 %>%
  #   mutate(obs_time_failureCR = pmin(Time_Failure1, Time_Failure2),
  #          obs_time = pmin(Time_Failure1,Time_Failure2,Time_Censor,Time_Tau),
  #          status = ifelse(obs_time == Time_Failure1, 1,
  #                 ifelse(obs_time == Time_Failure2, 2,
  #                        0))
  #          )
  # # # Calculate obs_time as the minimum of failure times and censoring time
  # # df_times$obs_time <- do.call(pmin, c(df_times[c("Time_Failure1",
  # #                                         "Time_Failure2",
  # #                                         "Time_Censor",
  # #                                         "Time_Tau")], na.rm = TRUE))
  # #
  # # # Create status variable based on obs_time and different failure times
  # # df_times$status <- ifelse(df_times$obs_time == df_times$Time_Failure1, 1,
  # #                       ifelse(df_times$obs_time == df_times$Time_Failure2, 2,
  # #                              0))
  # #
  # # View(df_times)
  # df_surv1 = df_times %>%
  #   # overall death indicator (failure from any cause = 1; 0 o/w)
  #   mutate(indD = ifelse(obs_time_failureCR <= pmin(Time_Censor, Time_Tau), 1, 0))
  # df_surv1$obs_time = df_surv1$obs_time*365
  # df_surv2 = inner_join(df_cov, df_surv1, by = "ID")
  # print(head(df_surv1))
  # print(head(df_cov))
  # print(head(df_surv2))
  # message("censoring rate: ", sum(df_surv2$status==0)/N,
  #         "\nfailure1 rate: ", sum(df_surv2$status==1)/N,
  #         "\nfailure2 rate: ", sum(df_surv2$status==2)/N,
  #         "\noverall failure rate: ", sum(df_surv2$status!=0)/N)
  # print(quantile(failure_t1))
  # print(quantile(failure_t2))
  # plot(failure_t1, main = "failure from cause 1")
  # plot(failure_t2, main = sprintf("failure from cause 2 for sim#%s", sim))
  # return(df_surv2)
  return(df_multi)
}

gdata_RE <- function(N=10,
                     G=2,
                  lambda_0D=0.1,lambda_0R=4,beta_R=log(2),beta_D=log(3),
                  omega_D=0.2,omega_R=0.1, gamma_D = 0.2, gamma_R = 0.1,
                  ztype=0, zparam=0.5,ctype=1,cparam=2,
                  gaptype=0,gapparam1=0.2,gapparam2=0.25,#rho;
                  num_A=2,tau=10, seed1=2023){

  # ztype indicates the distribution for covariate z: 0=normal(0,1),1=binary(zparam),2=uniform(0,1)
  # ctype indicates the distribution for censoring (0=exponential,1=uniform,
  #                                                 9=no censoring)
  set.seed(seed1)
  name = sprintf("Dataset_N%s_G%s_A%s_lambda0D%sBetaD%somegaD%sgammaR%s_lambda0R%sBetaR%somegaR%sgammaR%s_rho1%s_rho2%s_tau%s",
                 N, G, num_A,
                 round(lambda_0D,1), round(beta_D,1), round(omega_D,1), round(gamma_D,1), round(lambda_0R,1), round(beta_R,1), round(omega_R,1), round(gamma_R,1),
                 gapparam1, gapparam2, tau)

  # generating covariates
  if (ztype == 0){z <- rnorm(N)} #print("Covariates follow N(0,1)");
  if (ztype == 1){z <- rbinom(N, 1, zparam)}
  if (ztype == 2){z <- runif(N)}
  # print(sprintf("covariates z: %s",z))

  # generating treatment
  A = rbinom(N, num_A-1, 0.5)
  # print(sprintf("treatments A: %s",A))
  df_cov = data.frame(ID = c(1:N), Cov = z, Trt = A)

  # generating censoring time
  if (ctype == 0){cc <- rexp(N,cparam)}
  if (ctype == 1){cc <- runif(N,min=0,max=tau)}
  # print(sprintf("censoring time cc: %s",cc))

  # generating gap time 1 using Gumbel bivariate exponential model
  #if alpha=0 then independent and don't need to do this.
  alpha1 = 4*gapparam1
  alpha2 = 4*gapparam2
  lambda_D = lambda_0D*exp(t(beta_D)%*%z + omega_D * A + A * t(gamma_D) %*% z) #no intercept for z b/c survival
  lambda_R = lambda_0R*exp(t(beta_R)%*%z + omega_R * A + A * t(gamma_R) %*% z)
  # View(lambda_D)
  # View(lambda_R)
  if (alpha1 == 0){
    print("Independence")
    tt = rexp(N,lambda_0D*exp(t(beta_D)*z))
  } else{
    print("============== Gumbel Bivariate Exponential ==============")
    # generate first gap time using marginal exp(lamR) distribution
    gaptime1 = rexp(N,lambda_R)
    # print(sprintf("gaptime1: %s", gaptime1))
    # print("Generating Failure Time")
    # generate conditional failure time (given gaptime1)
    tt <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha1,y_type=1,y=gaptime1)$tt)
    # print("Checking Plot")
    plot_check = Check_GumbelBiExp(N=N,tt=tt);
    # print(plot_check)
    # Generate G conditional gaptime values for each subject
    if (G>1){
      gaptimes <- generate_gaptime(N, lambda_D, lambda_R, alpha2, G);#print(gaptimes)
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
      mutate(Gap = as.numeric(gsub("gaptime", "", Gap))); #print(head(gap_df_longitudinal))
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
    # gaptime2 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime1)$tt)
    # gaptime3 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime2)$tt)
    # gaptime4 <- as.numeric(GumbelBiExp(N=N,lambda_D=lambda_D,lambda_R=lambda_R,alpha=alpha2,y_type=2,y=gaptime3)$tt)
    # gap_df = data.frame(ID = c(1:N), gaptime1 = gaptime1, gaptime2 = gaptime2, gaptime3 = gaptime3, gaptime4 = gaptime4)
  }

  #steps to put together in one dataset
  df_times = data.frame(ID = c(1:N), Time_Failure=tt, Time_Censor=cc, Time_Tau=tau)

  # survival dataset (1row/person)
  df_surv1 = df_times %>% mutate(obs_time = pmin(Time_Failure, Time_Censor, Time_Tau),
                                 indD = ifelse(Time_Failure <= pmin(Time_Censor, Time_Tau), 1, 0)
                                 )
  df_surv2 = inner_join(df_cov, df_surv1, by = "ID") %>%
    dplyr::select(ID, Cov, Trt, obs_time, indD)

  # recurrent dataset
  df_times_long = df_times %>%
    pivot_longer(cols = starts_with("Time_"), names_to = "Label", values_to = "Time") %>%
    mutate(Label = gsub("Time_", "", Label)) %>%
    as.data.frame();
  df_long = rbind(recurrent_df_longitudinal, df_times_long) %>%
    dplyr::select(ID, Label, Time) %>%
    group_by(ID) %>%
    arrange(ID, Time)
  # find min(failure, censoring, tau)
  df_min = df_long %>%
    group_by(ID) %>%
    mutate(failure_status = ifelse(Label == "Failure", 1, ifelse(Label == "Censor" | Label == "Tau", 0, 99))) %>%
    filter(Label %in% c("Censor", "Failure", "Tau")) %>%
    summarise(min = min(Time),
              failure_status_raw = failure_status[which.min(Time)],
              Label = Label[which.min(Time)])
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
    mutate(L_open = ifelse(FirstRow ==1, 0, lag(Time)), #open paranthases ("(")
           R_closed = Time, #closed brack ("]")
           IndR = ifelse(contains_recurrent == TRUE, 1, 0),
           IndD = ifelse(contains_failure == TRUE, 1, 0)) %>%
    rename(R_Label = Label.x) %>%
    dplyr::select(ID, L_open, R_closed, IndR, IndD, R_Label) %>%
    as.data.frame()
  #add in covariates
  dataset = merge(dataset0,df_cov, by = "ID") %>%
    dplyr::select(ID, L_open, R_closed, IndR, IndD, Cov, Trt, R_Label)
  print(name)
  list(dataset_recurrent=dataset,dataset_survival=df_surv2,name=name)
}

# function for generating failure time from Gumbel bivariate exponential distribution (1960)
GumbelBiExp <- function(N,lambda_D,lambda_R,alpha,y_type=1,y=y) {
  # y_type indicates what type of conditional y we want: 1 = survival time; 2 = gap time
  # if y_type=2 then we want to input gaptime as the y aka gaptime2|gaptime1.
  if (y_type==1){
    title = "Generated Failure Time (conditional on first gap time)"
  } else if (y_type==2){
    title = "Generated Gap Time (conditional on previous gap time)"
  } else{
    stop("Invalid y_type")
  }
  u = runif(N); #print(sprintf("u: %s", u))
  c = 1-2*exp(-lambda_R*y); #print(sprintf("1-2exp(-y): %s", round(c,3)))
  rootsqd = (u-1)/(c*alpha) + ((1+c*alpha)/(2*c*alpha))^2; #print(sprintf("root^2: %s",rootsqd))
  root = sqrt(rootsqd)
  expnegx_plus = (1+c*alpha)/(2*c*alpha) + root;
  expnegx_minus = (1+c*alpha)/(2*c*alpha) - root;
  # print(sprintf("Indicator that exp(-x) (with PLUS root) > 1: %s", (expnegx_plus > 1)))
  # print(sprintf("Indicator that exp(-x) (with PLUS root) <= 1: %s", (expnegx_plus <= 1)))
  expnegx = expnegx_minus*(expnegx_plus>1)+expnegx_plus*(expnegx_plus<=1)
  # print(sprintf("exp(-x): %s", round(expnegx,3)))
  tt = -log(expnegx)
  # print(sprintf("%s (N = %s) is: %s", title, N, round(tt,2)))
  list(tt=tt)
}

generate_gaptime <- function(N, lambda_D, lambda_R, alpha, G) {
  gaptime <- matrix(0, nrow = N, ncol = G)
  for (i in 2:G) { # starting with gaptime2
    # print(sprintf("Generating Gap Time %s", i))
    gaptime[,i] <- as.numeric(GumbelBiExp(N = N, lambda_D = lambda_D, lambda_R = lambda_R,
                                          alpha = alpha, y_type = 2, y = gaptime[i - 1])$tt)
    # print(sprintf("End of Gap Time %s", i ))
  }
  return(gaptime)
}

Check_GumbelBiExp <- function(N,tt){
  # print("Plotting Estimated Cumulative Hazard over Time from Gumbel")
  # print("sdfljsdlfjsldfjlsd")
  # print(tt)
  # print(rep(1,N))
  s1 <- survfit(Surv(as.numeric(tt), rep(1,N)) ~ 1)
  # print(s1)
  LAM1 = -log(s1$surv)[1:N-1]
  # print(LAM1)
  time1 = s1$time[1:N-1]
  # print(time1)
  df = data.frame(time1,LAM1)
  # print(head(df))
  # print(length(time1)); print(length(LAM1))
  plot1 = ggplot(data = df, aes(x = time1, y = LAM1)) +
    geom_line() +
    xlab("Time") +
    ylab("Estimated Cumulative Hazard") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggtitle("Estimated cumulative hazard over time: generated failure time") +
    theme_bw()
  return(plot1)
}

# time calculator
tt <- function(s, reset = FALSE, units = "auto"){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    result <- data.frame(begin = time.tmp, end = Sys.time(), elapsed = difftime(Sys.time(), time.tmp, units = units))
    if (reset) time.tmp <<- Sys.time()
    return(result)
  }
}


flowchart <- function(outputData) {
  n.stages = dim(outputData)[3]
  n.sample = dim(outputData)[1]
  result <- data.frame(stage = c(1:n.stages, "Total", "Percent"),
                       total = rep(NA, n.stages + 2), percentage = NA,
                       censored = NA, died = NA, nextStage = NA)
  for (stage in 1:n.stages) {
    result[stage, "total"]  = dim(outputData)[1]
    result[stage, "percentage"]  = round(result[stage, "total"] / n.sample,2)
    cens <- outputData[, "delta", stage] == 0
    # drop those censored
    outputData <- outputData[!cens,,, drop = FALSE]
    died <- outputData[, "gamma", stage] == 1
    # drop those died
    outputData <- outputData[!died,,, drop = FALSE]
    result[stage, "censored"]  = sum(cens, na.rm = TRUE)
    result[stage, "died"]      = sum(died, na.rm = TRUE)
    result[stage, "nextStage"] = sum(!died, na.rm = TRUE)
  }
  result[n.stages + 1, c("censored", "died")] <- apply(result[, c("censored", "died")], 2, sum, na.rm = TRUE)
  result[n.stages + 2, c("censored", "died")] <- round(result[n.stages + 1, c("censored", "died")]/n.sample, 2)
  # result[c("Total", "Percent"), 1:2]
  result
}


output2observable <- function(output, cumulative = FALSE, evaluate = FALSE) {
  stage = n.stages = 1
  n = dim(output)[1]
  nm = dimnames(output)[[2]]
  nm.covar = c(grep("Z[0-9]+", nm, value = TRUE))
  if (evaluate == TRUE){
    nm.stage = c("OS_eval", "CIF_eval", "action", nm.covar)
  } else{
    nm.stage =  c("event.time", "status", "action", nm.covar)
  }
  df <- data.frame(subject.id = output[, "subj.id"])
  df.i <- data.frame(output[, nm.stage])
  if (evaluate == TRUE){
    names(df.i) = c("OS_eval", "CIF_eval", "A", nm.covar)
  } else{
    names(df.i) =  c("event.time", "status", "A", nm.covar)
  }
  df <- cbind(df, df.i)
  df
}
output2observable_phase <- function(output, phase = NULL, cumulative = FALSE) {
  if (length(dim(output)) > 2){
    n.phases = dim(output)[3]
  } else{
    print("one phase")
    n.phases = 1
  }

  if (is.null(phase)) {
    phase = 1:n.phases
  } else if (cumulative) {
    phase = 1:max(phase)
  }

  print(dim(output)[1])
  n = dim(output)[1]
  print(n)
  nm = dimnames(output)[[2]]
  nm.covar = c(grep("Z[0-9]+", nm, value = TRUE))
  nm.phase =  c("event.time", "delta", "action", nm.covar)

  if (length(dim(output)) > 2){
    df <- data.frame(subject.id = output[, "subj.id", 1])
    for (i in phase) {
      if (i == 1){ # phase 1
        df.i <- data.frame(output[, nm.phase, i])
        names(df.i) <- c(paste(c("T", "D", "A"), i, sep = "."), nm.covar)
      } else{
        df.i <- data.frame(output[, c("event.time", "delta", "action"), i])
        names(df.i) <- paste(c("T", "D", "A"), i, sep = ".")
      }
      df <- cbind(df, df.i)
    }
  } else{
    print("one phase")
    df <- data.frame(subject.id = output[, "subj.id"])
    i = 1
    print(nm.phase)
    print(head(output))
    df.i <- data.frame(output[, nm.phase])
    names(df.i) <- paste(c("T", "delta", "A", nm.covar), i, sep = "_")
    df <- cbind(df, df.i)
  }

  df$delta =
    df %>% dplyr::select(starts_with("delta_")) %>% as.matrix %>%
    apply(1, function(s) 1 - any(s == 0, na.rm = TRUE))
  df
}


# for (i in seq_along(rda_methods)) {
#   assign(paste0("skip.", rda_methods[i]), skip_method[i])
# }
# creating skip functions
assign_skip_function = function(methods_vector, skip_vector){
  for (i in seq_along(methods_vector)) {
    assign(paste0("skip.", methods_vector[i]), skip_vector[i], envir = .GlobalEnv)
  }
  return(skip_vector)
}

# CR02.Simulation_Summary.R function to select "XX_survival" or "XX_endpoint" for result.comb
# Generate the column names for each method
select_method_endpoints = function(method_vec, endpoint){
  method_variables <- lapply(method_vec,
                             function(method) {
                               paste0(method, "_", endpoint)
                             }) %>% unlist()
  return(method_variables)
}


# supplmentary functions
# plotting_indiv_curves_after_prediction.R in Scratch folder
indiv_predicted_Opt_curves = function(person, time_points, time_points_full,
                                      big_dataset_opt, endpoint, tau = 0){
  print(tau)
  if (tau == 1){
    tp = time_points_full
  }
  else {
    tp = time_points
  }
  opt_data = t(big_dataset_opt)[1:length(tp), person] %>%
    as.data.frame()
  rownames(opt_data) = NULL
  opt_data$tp = tp

  opt_plot = ggplot(opt_data, aes(x = tp, y = .)) +
    geom_step(color = 'violet', size = 1, linetype = 'solid') +
    # geom_point(color = 'blue', size = 2) +
    theme_minimal()
  if (tau == 1){
    last_tp = time_points[length(time_points)]
    opt_plot = opt_plot +
      geom_vline(xintercept = last_tp, color = 'red', linetype = 'solid') +
      labs(x = 'Time Points',
           y = sprintf("%s",endpoint),
           title = sprintf('Truncated %s for One Person Over Time', endpoint))
  } else{
    opt_plot = opt_plot +
      labs(x = 'Time Points',
           y = sprintf("%s",endpoint),
           title = sprintf('%s for One Person Over Time', endpoint))
  }
  return(opt_plot)
}

indiv_predicted_Trt_curves = function(person,
                                      time_points,
                                      time_points_full,
                                      dataset_trt0,
                                      dataset_trt1,
                                      endpoint,
                                      tau = 0){
  print(tau)

  if (tau == 1){
    tp = time_points_full
    last_tp = time_points[length(time_points)]
  } else{
    tp = time_points
  }
  trt0_data = as.data.frame(dataset_trt0)[1:length(tp), person]
  trt1_data = as.data.frame(dataset_trt1)[1:length(tp), person]
  trt_data = cbind(trt0 = trt0_data, trt1 = trt1_data) %>% as.data.frame()

  if (tau == 1){
    trt_plot = ggplot(trt_data, aes(x = tp)) +
      geom_vline(xintercept = last_tp, color = 'red', linetype = 'solid') +
      labs(x = 'Time',
           y = sprintf("%s", endpoint),
           title = sprintf('Truncated %s Function for ID=%s comparing Trt -1 vs Trt 1', endpoint, person))
    } else{
    trt_plot = ggplot(trt_data, aes(x = tp)) +
      labs(x = 'Time',
           y = sprintf("%s", endpoint),
           title = sprintf('%s Function for ID=%s comparing Trt -1 vs Trt 1', endpoint, person))
    }
  trt_plot = trt_plot +
    geom_step(aes(y = trt0, color = 'Trt -1'), size = 0.6, linetype = 'solid') +
    geom_step(aes(y = trt1, color = 'Trt 1'), size = 0.6, linetype = 'solid') +
    theme_minimal() +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_color_manual(name = "Treatment",
                       values = c("Trt -1" = "purple", #match the name of the element to the name in the aes
                                  "Trt 1" = "light blue")) #+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(trt_plot)
}

plotting_indiv_predicted_curves = function(person, time_points,time_points_full,
                                           big_dataset_opt_St, dataset_trt0_St, dataset_trt1_St,
                                           big_dataset_opt_CIF, dataset_trt0_CIF, dataset_trt1_CIF){
  # optimal Y for truncated time
  St_opt_plot = indiv_predicted_Opt_curves(person, time_points, time_points_full, big_dataset_opt_St, endpoint = "St")
  CIF_opt_plot = indiv_predicted_Opt_curves(person, time_points, time_points_full, big_dataset_opt_CIF, endpoint = "Cause1-CIF")

  # predicted by treatment for truncated time
  St_trt_plot = indiv_predicted_Trt_curves(person, time_points, time_points_full, dataset_trt0_St, dataset_trt1_St, endpoint = "St")
  CIF_trt_plot = indiv_predicted_Trt_curves(person, time_points, time_points_full, dataset_trt0_CIF, dataset_trt1_CIF, endpoint = "Cause1-CIF")

  # optimal Y for FULL time (to tau)
  St_opt_plot_tau = indiv_predicted_Opt_curves(person, time_points, time_points_full, big_dataset_opt_St, endpoint = "St", tau = 1)
  CIF_opt_plot_tau = indiv_predicted_Opt_curves(person, time_points, time_points_full, big_dataset_opt_CIF, endpoint = "Cause1-CIF", tau = 1)

  # predicted by treatment FULL time (to tau)
  St_trt_plot_tau = indiv_predicted_Trt_curves(person, time_points, time_points_full, dataset_trt0_St, dataset_trt1_St, endpoint = "St", tau = 1)
  CIF_trt_plot_tau = indiv_predicted_Trt_curves(person, time_points, time_points_full, dataset_trt0_CIF, dataset_trt1_CIF, endpoint = "Cause1-CIF", tau = 1)

  return(list(person = person,
              time_points = time_points,
              St_opt_plot = St_opt_plot,
              St_trt_plot = St_trt_plot,
              CIF_opt_plot = CIF_opt_plot,
              CIF_trt_plot = CIF_trt_plot,
              St_opt_plot_tau = St_opt_plot_tau,
              St_trt_plot_tau = St_trt_plot_tau,
              CIF_opt_plot_tau = CIF_opt_plot_tau,
              CIF_trt_plot_tau = CIF_trt_plot_tau))
}
message("End of F01.Simulation_Functions.R")


# End of script -------------------------------------------------------------

