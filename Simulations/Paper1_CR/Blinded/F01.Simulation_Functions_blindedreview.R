# Title: Libraries and Functions for Simulation Scripts
# Date: 07.02.2023

source("F00.Simulation_Libraries.R")
source("F02.ComparatorMethod_Functions.R")
source("F01.DynamicsCR.R")
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
  
  # message("ztype is:", ztype)
  if (ztype == 0) {
    # Generate random binary and continuous covariates for each variable
    # Randomly assign each covariate to binary or continuous with 50% chance
    # set.seed(zseed)
    cov_type <- sample(c(0, 1), ncov, replace = TRUE)  # 0 = continuous, 1 = binary
    # message("cov_type is:", cov_type)
    z <- matrix(0, nrow = N, ncol = ncov)   # Initialize covariate matrix 
    for (i in 1:ncov) {
      if (cov_type[i] == 1) {
        # Binary covariate (using binomial distribution)
        z[, i] <- rbinom(N, 1, 0.5)  # 50% probability for binary values
      } else {
        # Continuous covariate (normally distributed)
        z[, i] <- rnorm(N, mean = 0, sd = zparam)
      }
    }
  }
  if (ztype == 1) {
    # Generate binary covariates using binomial distribution
    z <- matrix(rbinom(N * ncov, 1, 0.5), nrow = N, ncol = ncov)  # 50% probability for binary values
  }
  if (ztype == 2) {
    # Generate continuous covariates using normal distribution
    z <- matrix(rnorm(N * ncov, mean = 0, sd = zparam), nrow = N, ncol = ncov)
  }
  if (ztype == 3) {
    # Generate continuous covariates using uniform distribution
    z <- matrix(runif(N * ncov), nrow = N, ncol = ncov)
  }
  if (is.null(colnames(z))) {
    colnames(z) <- paste0("Z", 1:ncov)
  }
  zzcov <<- z
  # View(zzcov)
  
  # generating censoring time
  if (ctype == 0){
    message("censoring: exp")
    # message("cparam:", cparam)
    cc <- rexp(N,cparam)
    # View(cc)
    # message("tau:", tau)
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
  df_multi <<- Dynamics(N=N, u1 = u1, u2=u2, u3 = u3,
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

  return(df_multi)
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

