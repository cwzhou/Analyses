# Title: Libraries and Functions for Simulation Scripts
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023

# Functions -----------------------------------------------------------------

gdata <- function(N=10,G=2,
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
  if (ztype == 0){print("Covariates follow N(0,1)"); z <- rnorm(N)} 
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
  View(lambda_D)
  View(lambda_R)
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
    print(plot_check)
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
  list(dataset_recurrent=dataset,name=name)
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


  


# End of script -------------------------------------------------------------
