num_cov = 1
# arg2 propensity
propensity <-   # (int), covariate (1~5)
  list(obs = list(beta.propensity = function(p) c(0, 1, -rep(0.5, p))),  # no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(0, 0, rep(0, p))))    # RCT

# arg1 beta
betas <- list(
  beta1 =
    list(beta.hazard0 = c(0, rep(1,num_cov)),              # (int), covariate (1)
         beta.hazard1 = c(0, c(2)),              # (int), covariate (1)
         beta.censor0 = c(-3, 0.2 * rep(1,num_cov)),      # (int), covariate (1)
         beta.censor1 = c(-3, 0.2 * (-1)))#,     # (int), covariate (1)
  # beta2 =  # moderate censoring rate
  #   list(beta.hazard0 = c(0, c(1,1,1, 1, 1)),              # (int), covariate (1~5)
  #        beta.hazard1 = c(0, c(2,2,1,-1,-1)),              # (int), covariate (1~5)
  #        beta.censor0 = c(-2, 0.2 * c(1,1,1, 1, 1)),      # (int), covariate (1~5)
  #        beta.censor1 = c(-2, 0.2 * c(1,1,1,-1,-1))),     # (int), covariate (1~5)
  # beta3 =  # common.beta1 + smaller p
  #   list(beta.hazard0 = c(0, -1, 1, 1),                        # (int), covariate (1~2)
  #        beta.hazard1 = c(0, -2, 2, 0),                        # (int), covariate (1~2)
  #        beta.censor0 = c(-3, 0.1, 0.2, 0.2),                  # (int), covariate (1~2)
  #        beta.censor1 = c(-3, 0.2, 0.2, -0.2)),                 # (int), covariate (1~2)
  # beta4 =  # common.beta1 + larger p
  #   list(beta.hazard0 = c(0, rep(1,8),  rep(0,2)),    # (int), covariate (1~10)
  #        beta.hazard1 = c(0, rep(2,4), rep(-1,4), rep(-2,2)),    # (int), covariate (1~10)
  #        beta.censor0 = c(-3, 0.2 * c(rep(1,4), rep(1,3), rep(0,3))),   # (int), covariate (1~10)
  #        beta.censor1 = c(-3, 0.2 * c(rep(1,4), rep(-1,3), rep(0,3))))   # (int), covariate (1~10)
)

p.list <- lapply(betas, function(x) length(x$beta.hazard0) - 2)
# setting1 n.boot 50 / n 300
setting = c(arg = list(arg), default, betas[[arg1]], p = p.list[[arg1]],
            propensity[[arg2]], size[[arg3]], crit[[2]])
# setting = c(arg = list(arg), default, betas[[arg1]], p = p.list[[arg1]],
#             propensity[[arg2]], size[[arg3]], crit[[arg4]])


gdata_CR <- function(N=10,
                     predHazardFn, predPropensityFn, # list of predictor functions
                     ztype=0, zparam=0.5, # covariates
                     ctype=1,cparam=2, # censoring
                     num_A=2, tau=10, seed1=2023){
  
  # ztype indicates the distribution for covariate z: 0=normal(0,1),1=binary(zparam),2=uniform(0,1)
  # ctype indicates the distribution for censoring (0=exponential,1=uniform,
  #                                                 9=no censoring)
  set.seed(seed1)
  name = sprintf("Dataset_N%s_G%s_A%s_lambda0D%sBetaD%somegaD%sgammaR%s_lambda0R%sBetaR%somegaR%sgammaR%s_rho1%s_rho2%s_tau%s",
                 N, G, num_A, 
                 round(lambda_0D,1), round(beta_D,1), round(omega_D,1), round(gamma_D,1), round(lambda_0R,1), round(beta_R,1), round(omega_R,1), round(gamma_R,1),
                 gapparam1, gapparam2, tau)
  
  # generating covariates
  if (ztype == 0){print("Covariates follow N(0,1)"); z <- rnorm(N)} 
  if (ztype == 1){sprintf("Covariates follow Bern(%s)", zparam); z <- rbinom(N, 1, zparam)}
  if (ztype == 2){print("Covariates follow Unif(0,1)"); z <- runif(N)}
  # print(sprintf("covariates z: %s",z))
  
  # generating treatment
  A = rbinom(N, num_A-1, 0.5)
  # print(sprintf("treatments A: %s",A))
  df_cov = data.frame(ID = c(1:N), Cov = z, Trt = A)
  print(head(df_cov))
  
  # generating censoring time
  if (ctype == 0){print("Censoring from exponential distribution");cc <- rexp(N,cparam)}
  if (ctype == 1){print("Censoring from Uniform distribution");cc <- runif(N,min=0,max=tau)}
  # print(sprintf("censoring time cc: %s",cc))
  
  # Competing Risks
  u2 = runif(N)
  
  
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