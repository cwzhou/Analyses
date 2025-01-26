#### This script is run automatically from CR02.Simulation_Run.R. DO not run alone.

# library(parallel)
setting_seed = 0

start_time = Sys.time()
column_names = NULL
realistic_train_truth_list = realistic_test_truth_list = list()

if (criterion_phase1 == "mean.prob.combo"){
  dtr_criterion = "surv.mean"
} else if (criterion_phase1 == "prob"){
  dtr_criterion = "surv.prob"
} else{
  dtr_criterion = "mean"
}

overall_survival_val.fn <- function(data){
  # mean truncated overall survival time (from any cause) or
  # mean overall survival probability at t0 = crit.eval.os
  mean(data$OS_eval, na.rm = TRUE)
}
endpoint_val.fn <- function(data) {
  # mean truncated years lost due to cause 1 (priorty cause) or
  # mean CIF probability at t0 = crit.eval.cif
  mean(data$CIF_eval, na.rm = TRUE)
}

if (criterion_phase1 %in% c("mean", "area", "surv.area")) {
  crit.eval.os = "mean"
  splitRule = "mean_surv"# for itrSurv Phase 1
  # overall_survival_val.fn <- function(data){
  #   mean(subset(data$event.time, data$status > 0),
  #        na.rm = TRUE)
  #   }
} else if (criterion_phase1 %in% c("prob", "mean.prob.combo")) {
  crit.eval.os = "prob"
  splitRule = "logrank_surv" # for itrSurv Phase 1
  # overall_survival_val.fn <- function(data) {
  #   mean(subset(data$event.time,
  #               data$status > 0) >= crit.value_phase1,
  #        na.rm = TRUE)}
}
# NOTE: below is assuming endpoint = CR aka we want mean LESS THAN critvalue
if (criterion_phase2 %in% c("mean", "area", "surv.area")) {
  crit.eval.cif = "mean"
} else if (criterion_phase2 %in% c("prob", "mean.prob.combo")) {
  crit.eval.cif = "prob"
  # endpoint_val.fn <- function(data, priority_cause) {
  #   mean(subset(data$event.time,
  #               # MAKE SURE THIS IS RIGHT since survival is >=
  #               data$status == priority_cause) <= crit.value_phase2,
  #        na.rm = TRUE)}
}

print(filename)
filename.tmp <- gsub("\\.rds", "_tmp.rds", filename)
if (local == 0){
  extracted_part <- gsub(".*/([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})/(.*).rds", "\\2", filename)
  filename.true2 <- sprintf("%s/%s_true2.rds", dir_rds_tmp, extracted_part)
} else{
  filename.true2 <- gsub("\\.rds", "_true2.rds", filename)
}

# Generate column names based on methods for result
column_names <- lapply(all_methods, function(method) {
  if (method == "obs" | method == "csk" | method == "pmcr" | method == "aipwe") {
    c(
      paste(method, "survival", sep = "_"),
      paste(method, "endpoint", sep = "_"),
      paste("time", method, sep = ".")
    )
  } else {
    c(
      paste(method, "survival", sep = "_"),
      paste(method, "endpoint", sep = "_"),
      paste(method, "n_phase2", sep = "_"),
      paste("time", method, sep = ".")
    )
  }
})
# Flatten the list of column names
column_names <- unlist(column_names)
# Define a custom sorting function
custom_sort <- function(names) {
  order(
    !grepl("_survival$", names),
    !grepl("_endpoint$", names),
    !grepl("_n_phase2$", names),
    !grepl("^time_", names),
    !grepl("^czmk_", names),
    !grepl("^csk_", names),
    !grepl("^pmcr_", names),
    !grepl("^aipwe_", names),
    !grepl("^zom_", names),
    !grepl("^obs_", names),
    names
  )
}
sorted_column_names <- column_names[custom_sort(column_names)]
# print(sorted_column_names)
# for r00
colnames1 = lapply(all_methods, function(method) {
  c(
    paste(method, "training_os_trtprop", sep = "_"),
    paste(method, "training_cause1_trtprop", sep = "_"),
    paste(method, "testing_os_trtprop", sep = "_")
  )
})
# Flatten the list of column names
column_names1 <- unlist(colnames1)

if (endpoint != "CR" | generate_failure_method != "fine_gray"){
  message("setting cause1prob to NULL b/c not in fine-gray setting")
  cause1prob = NULL
}
if (endpoint == "CR"){
  arg_list =   list(
    N = n, tau = tau, # structural parameters
    ztype = 2, #uniform covariates
    ctype=ctype,
    cparam=censor_rate,
    censor_min=censor_min,
    censor_max=censor_max,
    ncov = ncov,
    # M = M, #number of causes # NOTE WE CURRENTLY ONLY HAVE GDATA_CR SUPPORT 2 CAUSES b/c of Fine-Gray simulation setting
    mass_p = cause1prob,
    generate_failure_method = generate_failure_method,
    predHazardFn = predHazardFn,
    predPropensityFn = predPropensityFn, # list of predictor functions,
    priority_cause = priority_cause
  )
} else{
  list(
    N = n, tau = tau, tick = 0.01,  # structural parameters
    censor_rate = censor_rate,
    ncov = ncov, #M = M,
    predHazardFn = predHazardFn,
    predPropensityFn = predPropensityFn, # list of predictor functions
    hidden_data = TRUE, printFlag = FALSE
  )
}

nodesize = 5
mindeath = round(sqrt(c(nodesize)), 0)
tau = tau

arg.obs <- arg.czmk <- arg.csk <-
  arg.pmcr <- arg.aipwe <-
  arg.obs.no.censor <- arg.zom <- arg_list

# arg.czmk$priority_cause <- arg.csk$priority_cause <- arg.zom$priority_cause <- priority_cause

# for NOW only: no censoring
# arg.obs$ctype <- 99 # comment this out when we want training data to have censoring

# evaluation datasets dont have censoring
arg.czmk$ctype <- arg.csk$ctype <-
  arg.pmcr$ctype <- arg.aipwe$ctype <-
  arg.zom$ctype <- arg.obs.no.censor$ctype <- 99 # no censoring

arg.czmk$N <- arg.csk$N <-
  arg.pmcr$N <- arg.aipwe$N <-
  arg.obs.no.censor$N <- arg.zom$N <- n.eval

arg.czmk$crit.eval.cif <- arg.csk$crit.eval.cif <-
  arg.pmcr$crit.eval.cif <- arg.aipwe$crit.eval.cif <-
  arg.obs.no.censor$crit.eval.cif <- arg.zom$crit.eval.cif <- crit.eval.cif
arg.czmk$crit.eval.os <- arg.csk$crit.eval.os <-
  arg.pmcr$crit.eval.os <- arg.aipwe$crit.eval.os <-
  arg.obs.no.censor$crit.eval.os <- arg.zom$crit.eval.os <- crit.eval.os

#tau/2 is because of itrSurv::VerifySurvivalTime.R (line 53 etc) used in criticalValue.R
if (is.null(crit.value_phase2) || is.na(crit.value_phase2) || length(crit.value_phase2) == 0){
  t0.cif = tau/2;
  t0_pmcr = 0.1 #tau/2
  t0_aipwe = 0.1 # these need to be tuned...
} else{
  t0.cif = crit.value_phase2;
  t0_pmcr = 0.1 #crit.value_phase2
  t0_aipwe = 0.1 # this needs to be tuned..
}
if (is.null(crit.value_phase1) || is.na(crit.value_phase1) || length(crit.value_phase2) == 0){
  t0.os = tau/2;
} else{
  t0.os = crit.value_phase1;
}
arg.czmk$t0.cif <- arg.csk$t0.cif <-
  arg.pmcr$t0.cif <- arg.aipwe$t0.cif <-
  arg.obs.no.censor$t0.cif <- arg.zom$t0.cif <- t0.cif;
arg.czmk$t0.os <- arg.csk$t0.os <-
  arg.pmcr$t0.os <-arg.aipwe$t0.os <-
  arg.obs.no.censor$t0.os <- arg.zom$t0.os <- t0.os;

arg.czmk$evaluate <- arg.csk$evaluate <-
  arg.pmcr$evaluate <- arg.aipwe$evaluate <-
  arg.obs.no.censor$evaluate <- arg.zom$evaluate <- TRUE

# Define the function that runs a single simulation
# run_simulation <- function(sim){
message("starting run_simulation")
### simulation
r00 = data.frame(sim = n.sim) # this is for parallelizaitn: data.frame(sim.no = sim) #this is for sim for loop: data.frame(sim = n.sim)
r00[column_names1] <- NA

# Create the data frame with the 'sim' column
result <- data.frame(sim = 1:n.sim) # this is for parallelizaitn: data.frame(sim.no = sim) #this is for sim for loop: data.frame(sim = 1:n.sim)
# Add columns to the result data frame
result[sorted_column_names] <- NA
attr(result, "criterion_phase1") <- list(criterion = criterion_phase1, crit.value = crit.value_phase1)
attr(result, "criterion_phase2") <- list(criterion = criterion_phase2, crit.value = crit.value_phase2)

trt_result <- data.frame(rep(1:n.sim, each = n.methods), # this is for parallelization: sim.no = rep(sim, each = n.methods), #this is for sim for loop: rep(1:n.sim, each = n.methods),
                         method = NA,
                         surv_A = NA, surv_B = NA,
                         endpoint_A = NA, endpoint_B = NA) # each = 3 means repeating b/c theres 3 methods rn (obs, czmk, zom)

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

# flow: obs (1) -> optimal (n.mc), obs.no.censor (n.mc)
true2 = true3 = list()
print(Sys.time())
for (sim in 1:n.sim) {
  
  cat("\n\n#################################")
  cat("\n######### Simulation ",sim, "#########")
  cat("\nEndpoint: ",endpoint)
  if (endpoint == "CR"){
    cat("\nNumber of CR: ",M)
  }
  cat("\n#################################\n")
  
  # if (sim == 1){
  #   sim_count = sim
  #   prev_count = 0
  # } else{
  #   sim_count = sim*2+prev_count
  #   prev_count = prev_count + 1
  # }
  
  train_seed = sim*10000 + 123
  # if (setting_seed == 1){
  #   arg.obs$seed1 <- train_seed
  #   arg.czmk$seed1 <- arg.obs.no.censor$seed1 <-
  #   arg.zom$seed1 <- arg.csk$seed1 <- train_seed + 10 #eval seed (generate same covariates)
  # }
  
  cat ("1. Data for Simulation",sim, ":",generate_failure_method,"\n")
  tt(1)
  
  # obs.data
  set.seed(train_seed + 0)
  u1_train = runif(n)
  u1_test = runif(n.eval)
  u2_train = runif(n)
  u2_test = runif(n.eval)
  u3_train = runif(n)
  u3_test = runif(n.eval)
  
  arg.obs$u1 = u1_train; arg.obs$u2 = u2_train; arg.obs$u3 = u3_train;
  arg.czmk$u1 <- arg.csk$u1 <- arg.pmcr$u1 <- arg.aipwe$u1 <- arg.zom$u1 <- arg.obs.no.censor$u1 <- u1_test
  arg.czmk$u2 <- arg.csk$u2 <- arg.pmcr$u2 <- arg.aipwe$u2 <- arg.zom$u2 <- arg.obs.no.censor$u2 <- u2_test
  arg.czmk$u3 <- arg.csk$u3 <- arg.pmcr$u3 <- arg.aipwe$u3 <- arg.zom$u3 <- arg.obs.no.censor$u3 <- u3_test
  # print(head(u1_train,1));print(head(u1_test,1));print(head(u2_train,1));print(head(u2_train,1))
  arg.test0 = arg.obs.no.censor;
  arg.test0$policy = "test_action0"
  arg.test1 = arg.obs.no.censor;
  arg.test1$policy = "test_action1"
  arg.train0 = arg.obs;
  arg.train0$policy = "test_action0"
  arg.train1 = arg.obs;
  arg.train1$policy = "test_action1"
  arg.train0$crit.eval.cif <- arg.train1$crit.eval.cif <- crit.eval.cif
  arg.train0$t0.cif <- arg.train1$t0.cif <- crit.value_phase1
  arg.train0$crit.eval.os <- arg.train1$crit.eval.os <- crit.eval.os
  arg.train0$t0.os <- arg.train1$t0.os <- crit.value_phase2
  
  message("using train_seed for generating training data")
  set.seed(train_seed)
  obs.data <- do.call(gdata_CR, arg.obs);
  obs_1 <<- obs.data
  
  tmp_name = sprintf("meantimes_sim%s", sim)
  # # we want to have higher mean for action = 1
  overall = obs_1 %>%
    group_by(action) %>%
    summarize(mean = mean(event.time),
              n=n())
  # # we want to have higehr mean for action = -1 for status = 1
  cause = obs_1 %>%
    group_by(status,action) %>%
    summarize(mean = mean(event.time),
              n=n())
  m = list(overall = overall, cause = cause)
  assign(tmp_name, m)
  
  arg.obs_tmp = arg.obs;arg.obs_tmp$full_data = 1
  obs_1_tmp = do.call(gdata_CR, arg.obs_tmp)
  
  # View(obs.data)
  # print(head(obs_1$event.time))
  message("data.df")
  data.df <<- obs.data %>%
    mutate(D.0 = ifelse(status> 0,1,0), # event from any cause indicator
           D.1 = ifelse(status==1,1,0), # event from cause 1 indicator
           D.2 = ifelse(status==2,1,0), # event from cause 2 indicator
           obs_time = event.time,
           Trt = action) %>%
    dplyr::select(obs_time, Trt, status, D.0, D.1, D.2, contains("Z")) %>%
    arrange(obs_time) # MUST SORT ASCENDING ORDER FOR OBS_TIME
  # dplyr::select(-c(Time_Censor, Time_Failure1, Time_Failure2, Time_Tau, obs_time_failureCR,indD))
  
  # create unique observed failure times for survival phase 1
  timePointsSurvival = data.df %>%
    filter(D.0 == 1) %>%
    dplyr::select(obs_time) %>%
    unlist(use.names = FALSE) %>%
    unique()
  # for CR only
  timePointsEndpoint = timePointsSurvival # we can do this bc CR is subset of OS failure times
  # create unique observed endpoint times for endpoint phase 2
  
  data.df.pmcr <<- data.df %>%
    dplyr::select(obs_time, Trt, status, starts_with("Z")) %>%
    mutate(Trt = ifelse(Trt == -1, 0, Trt)) # Trt has to be 0/1
  
  data.df.aipwe <<- aipwe_data_format(data.df)
  
  result[sim, "training_percent.censor"] <- mean(data.df$status==0, na.rm = TRUE) #result["training_percent.censor"]
  # result[paste0("training_cause.", 1:n.causes)]
  result[sim, paste0("training_cause.", 1:n.causes)] <-
    sapply(1:n.causes, function(s) mean(data.df$status == s))
  
  # # # obs policy value
  message("using train_seed+10 to generate obs testing data")
  set.seed(train_seed + 10)
  obs.data.rep <- do.call(gdata_CR, arg.obs.no.censor) # no censoring for eval sets
  rep_obs <<- obs.data.rep
  
  # result["obs_survival"]
  # result["obs_endpoint"]
  result[sim, "obs_survival"] <- overall_survival_val.fn(obs.data.rep)#val.fn_phase1(obs.data.rep$event.time)
  result[sim, "obs_endpoint"] <- endpoint_val.fn(obs.data.rep)
  
  # # A = -1; B = 1
  # # within method, mean survival time from those who have action -1 and mean survival time for those who have action 1
  # A_mean = obs.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
  # B_mean = obs.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
  # A_mean_cr = obs.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
  # B_mean_cr = obs.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
  # # AB_diff = A_mean - B_mean
  # trt_result[sim_count, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("obs", A_mean, B_mean, A_mean_cr, B_mean_cr)
  # rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
  
  result[sim, "time.obs"] <- tt(2, reset = TRUE, unit = "min")["elapsed"] #result["time.obs"]
  # # print(flowchart(obs.data$output)) # if i fix this in 02_Simulation_Functions.R for CR then I can use, otherwise delete.
  # # transforming data from an array format to a data.frame format
  # # data.df = output2observable(obs.data)
  rm(obs.data.rep); gc()
  
  set.seed(train_seed + 10)
  test0.data.rep <<- do.call(gdata_CR, arg.test0)
  set.seed(train_seed + 10)
  test1.data.rep <<- do.call(gdata_CR, arg.test1)
  set.seed(train_seed)
  train0.data <<- do.call(gdata_CR, arg.train0)
  set.seed(train_seed)
  train1.data <<- do.call(gdata_CR, arg.train1)
  
  cat("\n******************************\n")
  # estimation
  cat ("2. czmk for simulation", sim, "\n")
  if (!skip.czmk) {
    if ("package:dtrSurv" %in% search()) {
      detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
    }
    cat ("  2. czmk - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    # new package itrSurv
    priority_vector <- c(0, priority_cause)
    # Use lapply to create a list of sprintf statements
    models_itr <- lapply(priority_vector, function(x) {
      paste0(sprintf("Surv(obs_time, D.%s) ~ ", ifelse(x == 0, "0", as.character(x))),
             paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
    })
    # print(crit.value_phase1)
    arg.czmk2 = list(data = data.df,
                     endPoint = "CR",
                     yName = "obs_time",
                     txName = paste("Trt"),
                     epName = "status", # status indicator (0 = censored, 1 = failure from cause 1 (PC), 2 = failure from cause 2)
                     models = models_itr,
                     timePointsSurvival = timePointsSurvival,
                     timePointsEndpoint = timePointsEndpoint,
                     tau = tau,
                     criticalValue1 = criterion_phase1,
                     criticalValue2 = criterion_phase2,
                     tol1 = tol1,
                     evalTime = crit.value_phase1,
                     splitRule1 = splitRule,
                     splitRule2 = "gray_cr",
                     ERT = TRUE,
                     uniformSplit = TRUE,
                     replace = FALSE,
                     randomSplit = 0.2,
                     nTree = 300,
                     pooled = FALSE,
                     stratifiedSplit = 0.1)
    set.seed(train_seed + 1)
    optimal.czmk <- do.call(itrSurv::itrSurv, c(arg.czmk2,
                                                list(mTry = sqrt(ncov),
                                                     nodeSize = nodesize,
                                                     minEvent = mindeath)))
    czmk.error <- class(optimal.czmk)[1] == "try-error"
    arg.czmk$policy <- if (!czmk.error) optimal.czmk
    # View(arg.czmk$policy)
    policy_czmk <<- arg.czmk$policy
    #training
    predd_surv_czmk <<- predd_surv
    predd_ep_czmk <<- predd_ep
    # print(arg.czmk$policy@phaseResults[["EndPointPhase2Results"]]@optimal@optimalTx)
    if (!czmk.error) result[sim, "czmk_n_phase2"] <- mean(arg.czmk$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) #result["czmk_n_phase2"]
    rm(optimal.czmk); gc()
    
    cat ("  \n 2. czmk - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_czmk <<- arg.czmk
    if (!czmk.error) {
      set.seed(train_seed + 10)
      czmk.data.rep <- do.call(gdata_CR, arg.czmk)
      predd_surv_czmk_eval <<- predd_surv
      predd_ep_czmk_eval <<- predd_ep
      rep_czmk <<- czmk.data.rep
      
      # if (criterion_phase1 == "mean.prob.combo"){
      #   # message("checking survival curve for first 3 timepoints for everyone")
      #   # View(as.data.frame(docalling1[["optimal"]]@optimalY[, 1:3]))
      #   indiv_St_opt <<- as.data.frame(docalling1[["optimal"]]@optimalY[,])
      #   indiv_St_trt0 <<- docalling1[["value"]][["Prob"]][[1]]
      #   indiv_St_trt1 <<- docalling1[["value"]][["Prob"]][[2]]
      #   # message("checking priority-cause CIF curve for first 3 timepoints for everyone")
      #   indiv_CIF_opt <<- as.data.frame(docalling2[["optimal"]]@optimalY[,])
      #   indiv_CIF_trt0 <<- docalling2[["value"]][["Prob"]][[1]]
      #   indiv_CIF_trt1 <<- docalling2[["value"]][["Prob"]][[2]]
      #   # plotting_indiv_curves_after_prediction.R
      # } else{
      # message("checking survival curve for first 3 timepoints for everyone")
      # View(as.data.frame(docalling1[["optimal"]]@optimalY[, 1:3]))
      indiv_St_opt <<- as.data.frame(docalling1[["optimal"]]@optimalY[,])
      indiv_St_trt0 <<- docalling1[["value"]][["Func"]][[1]]
      indiv_St_trt1 <<- docalling1[["value"]][["Func"]][[2]]
      # message("checking priority-cause CIF curve for first 3 timepoints for everyone")
      indiv_CIF_opt <<- as.data.frame(docalling2[["optimal"]]@optimalY[,])
      indiv_CIF_trt0 <<- docalling2[["value"]][["Func"]][[1]]
      indiv_CIF_trt1 <<- docalling2[["value"]][["Func"]][[2]]
      # plotting_indiv_curves_after_prediction.R
      # }
      
      
      
      result[sim, "czmk_survival"] = overall_survival_val.fn(czmk.data.rep) #result["czmk_survival"]
      result[sim, "czmk_endpoint"] = endpoint_val.fn(czmk.data.rep) #result["czmk_endpoint"]
      # result[sim, "czmk_percent.censor"] <- mean(czmk.data.rep$status==0, na.rm = TRUE)
      # result[sim, paste0("czmk_cause.", 1:n.causes)] <-
      #   sapply(1:n.causes, function(s) mean(czmk.data.rep$status == s))
    }
    
    # # within czmk method, mean survival time from those who have action -1 and mean survival time for those who have action 1
    # A_mean = czmk.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
    # B_mean = czmk.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
    # A_mean_cr = czmk.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
    # B_mean_cr = czmk.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
    # trt_result[sim_count+1, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("czmk", A_mean, B_mean, A_mean_cr, B_mean_cr)
    # rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
    
    result[sim, "time.czmk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.czmk"]
    arg.czmk$policy <- NULL; gc()
    rm(czmk.data.rep); gc()
  }
  
  cat("\n******************************\n")
  # estimation
  cat("3. csk - Cho et al for Simulation",sim, ":",generate_failure_method,"\n")
  n.stages = 1
  if (!skip.csk) {
    library(dtrSurv)
    cat ("  3. csk - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    # dtrSurv (Cho et al)
    models_dtr <- paste0("Surv(obs_time, D.0) ~ ",
                         paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
    arg.csk2 = list(data = data.df,
                    txName = paste("Trt"),
                    models = models_dtr,
                    usePrevTime = TRUE, tau = tau, timePoints = "uni", nTimes = 100,
                    criticalValue = dtr_criterion, evalTime = crit.value_phase1,
                    # splitRule = NULL,
                    splitRule = ifelse(dtr_criterion == "mean", "mean", "logrank"),
                    ERT = TRUE, uniformSplit = TRUE, replace = FALSE,
                    randomSplit = 0.2, nTree = 300,
                    mTry = rep(sqrt(ncov), n.stages),
                    pooled = FALSE,
                    stratifiedSplit = 0.1)
    set.seed(train_seed + 2)
    optimal.csk <- do.call(dtrSurv::dtrSurv, c(arg.csk2, list(nodeSize = nodesize, minEvent = mindeath )))
    csk.error <- class(optimal.csk)[1] == "try-error"
    arg.csk$policy <- if (!csk.error) optimal.csk
    policy_csk <<- arg.csk$policy
    rm(optimal.csk); gc()
    
    cat ("  \n 3. csk - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_csk <<- arg.csk
    # if (!csk.error) {csk.data.rep <- do.call(multiStageDynamics, arg.csk)}
    # if (!csk.error) result[sim, "csk"] <- val.fn(csk.data.rep$summary$cumulative.event.time)
    
    if (!csk.error) {
      set.seed(train_seed + 10)
      csk.data.rep <- do.call(gdata_CR, arg.csk)
      rep_csk <<- csk.data.rep
      result[sim, "csk_survival"] = overall_survival_val.fn(csk.data.rep) #result["csk_survival"]
      result[sim, "csk_endpoint"] = endpoint_val.fn(csk.data.rep) # result["csk_endpoint"]
      # result[sim, "csk_percent.censor"] <- mean(csk.data.rep$status==0, na.rm = TRUE)
      # result[sim, paste0("csk_cause.", 1:n.causes)] <-
      #   sapply(1:n.causes, function(s) mean(csk.data.rep$status == s))
    }
    # # within csk method, mean survival time from those who have action -1 and mean survival time for those who have action 1
    # A_mean = csk.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
    # B_mean = csk.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
    # A_mean_cr = csk.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
    # B_mean_cr = csk.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
    # # AB_diff = A_mean - B_mean
    # trt_result[sim_count+2, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("csk", A_mean, B_mean, A_mean_cr, B_mean_cr)
    # rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
    
    result[sim, "time.csk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.csk"]
    arg.csk$policy <- NULL; gc()
    rm(csk.data.rep); gc()
    if ("package:dtrSurv" %in% search()) {
      detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
    }
  }
  
  
  cat("\n******************************\n")
  # estimation
  cat ("4. pmcr for Simulation",sim, ":",generate_failure_method,"\n")
  if (!skip.pmcr) {
    cat ("  4. pmcr - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    models_pmcr <- paste0("Trt ~ ",
                          paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
    arg.pmcr2 = list(Time="obs_time",
                     Event="status",
                     formula=models_pmcr, # Trt has to be 0/1 for pmcr
                     data=data.df.pmcr, # b/c trt is 0/1
                     rgenoud=FALSE,
                     Restrict=FALSE, # required for unrestricted
                     propscore="logistic", # change as appropriate
                     t0=t0_pmcr) # t0 need to be tuned
    # Unrestricted regime
    message("unrestricted regime")
    set.seed(train_seed + 3)
    temp.unrestr.fit1 <- do.call(PMCR, arg.pmcr2)
    #range of alpha # alpha needs to be tuned
    alps<-c(temp.unrestr.fit1$Fbeta2[2],
            seq(round(temp.unrestr.fit1$Fbeta2[2],2)+0.01,
                round(temp.unrestr.fit1$Fbeta1[2],2)+0.03,
                0.01))
    # NEED TO DO: update t0 above to be appropriate as well
    #Restricted regime parameters
    arg.pmcr3 = arg.pmcr2;
    arg.pmcr3$Restrict = TRUE
    alp<-alps[2]; arg.pmcr3$alp = alp # alp needs to be tuned based on alps
    M_PMCR<-1e+5; arg.pmcr3$M = M_PMCR # pick a large M
    #Restricted regime
    message("restricted regime")
    set.seed(train_seed + 3)
    temp.restr.fit1 <- do.call(PMCR, arg.pmcr3)
    optimal.pmcr <- temp.restr.fit1$beta3 #optbetas for the linear decision rule
    
    pmcr.error <- class(optimal.pmcr)[1] == "try-error"
    arg.pmcr$policy = list()
    arg.pmcr$policy[[1]] = "pmcr"
    arg.pmcr$policy[[2]] <- if (!pmcr.error) optimal.pmcr
    arg.pmcr$policy[[3]] = t0_pmcr
    arg.pmcr$policy[[4]] = arg.pmcr3
    policy_pmcr <<- arg.pmcr$policy
    rm(optimal.pmcr); gc()
    
    cat ("  \n 4. pmcr - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_pmcr <<- arg.pmcr
    if (!pmcr.error) {
      set.seed(train_seed + 10)
      pmcr.data.rep <- do.call(gdata_CR, arg.pmcr)
      rep_pmcr <<- pmcr.data.rep
      result[sim, "pmcr_survival"] = overall_survival_val.fn(pmcr.data.rep) #result["pmcr_survival"]
      result[sim, "pmcr_endpoint"] = endpoint_val.fn(pmcr.data.rep) #result["pmcr_endpoint"]
      # result[sim, "pmcr_percent.censor"] <- mean(pmcr.data.rep$status==0, na.rm = TRUE)
      # result[sim, paste0("pmcr_cause.", 1:n.causes)] <-
      #   sapply(1:n.causes, function(s) mean(pmcr.data.rep$status == s))
    }
    # # within pmcr method, mean survival time from those who have action -1 and mean survival time for those who have action 1
    # A_mean = pmcr.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
    # B_mean = pmcr.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
    # A_mean_cr = pmcr.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
    # B_mean_cr = pmcr.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
    # trt_result[sim_count+3, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("pmcr", A_mean, B_mean, A_mean_cr, B_mean_cr)
    # rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
    
    result[sim, "time.pmcr"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.pmcr"]
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
  
  
  cat("\n******************************\n")
  # estimation
  cat ("4. aipwe for Simulation",sim, ":",generate_failure_method,"\n")
  if (!skip.aipwe) {
    cat ("  4. aipwe - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    # models_aipwe <- paste0("Trt ~ ",
    # paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
    arg.aipwe2 = list(data = data.df.aipwe, # b/c trt is 0/1
                      pp.v = ncov,#c(2,1,2),#c(1,1,1,2), #ncov,
                      tau1 = t0_aipwe, # t0 needs to be tuned
                      tune = c(0.001,0.01,0.5,1,
                               seq(0.1,400,
                                   length.out=16))
    )
    set.seed(train_seed + 4)
    aipwe.fit1 <- do.call(aipwe.fit, arg.aipwe2)
    optimal.aipwe <- aipwe.fit1 #eta0-eta_{ncov}
    aipwe.error <- class(optimal.aipwe)[1] == "try-error"
    arg.aipwe$policy = list()
    arg.aipwe$policy[[1]] = "aipwe"
    arg.aipwe$policy[[2]] <- if (!aipwe.error) optimal.aipwe
    arg.aipwe$policy[[3]] = t0_aipwe
    arg.aipwe$policy[[4]] = arg.aipwe2
    policy_aipwe <<- arg.aipwe$policy
    rm(optimal.aipwe); gc()
    
    cat ("  \n 4. aipwe - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_aipwe <<- arg.aipwe
    if (!aipwe.error) {
      set.seed(train_seed + 10)
      aipwe.data.rep <- do.call(gdata_CR, arg.aipwe)
      rep_aipwe <<- aipwe.data.rep
      result[sim, "aipwe_survival"] = overall_survival_val.fn(aipwe.data.rep) #result["aipwe_survival"]
      result[sim, "aipwe_endpoint"] = endpoint_val.fn(aipwe.data.rep) #result["aipwe_endpoint"]
    }
    
    # # within aipwe method, mean survival time from those who have action -1 and mean survival time for those who have action 1
    # A_mean = aipwe.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
    # B_mean = aipwe.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
    # A_mean_cr = aipwe.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
    # B_mean_cr = aipwe.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
    # trt_result[sim_count+3, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("aipwe", A_mean, B_mean, A_mean_cr, B_mean_cr)
    # rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
    
    result[sim, "time.aipwe"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.aipwe"]
    arg.aipwe$policy <- NULL; gc()
    rm(aipwe.data.rep); gc()
    # aipwe.fit(data_list = data.df.aipwe,
    #           pp.v = ncov,
    #           tau1 = 0.1,
    #           tune = tune)
  }
  
  #
  #   cat ("3. Goldberg & Kosorok - lm \n")
  #   if (!skip.gk) {
  #     cat ("  3. Goldberg & Kosorok - lm - Policy estimation \n")
  #     set.seed(train_seed + 2)
  #     optimal.gk.lm <-
  #       try(gk.separate(common.formula = formula.lm(ncov),
  #                       common.Tx.label = "A", stage.label = 1:n.stages, tau = tau,
  #                       data = data.df, stage.sep = "_", method = "lm", regress.prev.time = F))
  #     # optimal.gk.lm <- try(gk.Q(data = obs.data$output, est.fn = gk.lm, formula = formula.lm(ncov), tau = tau))
  #     gklm.error <- class(optimal.gk.lm)[1] == "try-error"
  #     arg.gk.lm$policy <- if (!gklm.error) optimal.gk.lm$survRF
  #     attr(arg.gk.lm$policy, "class") = "GKLM"
  #     rm(optimal.gk.lm); gc()
  #
  #     cat ("  3. Goldberg & Kosorok - lm - Evaluation \n")
  #     set.seed(train_seed + 10)
  #     if (!gklm.error) gklm.data.rep <- do.call(multiStageDynamics, arg.gk.lm)
  #     if (!gklm.error) result[sim, "gkLM"]    <- val.fn(gklm.data.rep$summary$cumulative.event.time)
  #     result[sim, "time.gkLM"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
  #     arg.gk.lm$policy <- NULL; gc()
  #     rm(gklm.data.rep); gc()
  #
  #     cat ("4. Estimation - Goldberg & Kosorok  - rf \n")
  #
  #     cat ("  4. Estimation - Goldberg & Kosorok  - rf - Policy estimation \n")
  #     set.seed(train_seed + 3)
  #     arg.gk.rf2 = list(common.formula = formula.rf(ncov),
  #                       common.Tx.label = "A", stage.label = 1:n.stages, tau = tau,
  #                       data = data.df, stage.sep = "_", regress.prev.time = F)
  #     optimal.gk.rf <-
  #       try(do.call(gk.separate, c(arg.gk.rf2, list(nodesize = nodesize, method = "rf"))))
  #     # optimal.gk.lm <- try(gk.Q(data = obs.data$output, est.fn = gk.lm, formula = formula.lm(ncov), tau = tau))
  #     gkrf.error <- class(optimal.gk.rf)[1] == "try-error"
  #     arg.gk.rf$policy <- if (!gkrf.error) optimal.gk.rf$survRF
  #     if (!gkrf.error) attr(arg.gk.rf$policy, "class") = "GKRF"
  #     rm(optimal.gk.rf); gc()
  #
  #     cat ("  4. Estimation - Goldberg & Kosorok  - rf - Evaluation \n")
  #     set.seed(train_seed + 10)
  #     if (!gkrf.error) gkrf.data.rep <- do.call(multiStageDynamics, arg.gk.rf)
  #     if (!gkrf.error) result[sim, "gkRF"]    <- val.fn(gkrf.data.rep$summary$cumulative.event.time)
  #     result[sim, "time.gkRF"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
  #     arg.gk.rf$policy <- NULL; gc()
  #     rm(gkrf.data.rep); gc()
  #   }
  #
  #   cat ("5. Estimation - Simoneau et al. \n")
  #   if (!skip.dw) {
  #     cat ("  5. Estimation - Simoneau et al. - Policy estimation \n")
  #     set.seed(train_seed + 4)
  #     optimal.dw <- try(dwSurv(data = obs.data$output, tau = tau))
  #     dw.error <- class(optimal.dw)[1] == "try-error"
  #     arg.dw$policy <- if (!dw.error) optimal.dw$dw.est
  #     rm(optimal.dw); gc()
  #
  #     cat ("  5. Estimation - Simoneau et al. - Evaluation \n")
  #     set.seed(train_seed + 10)
  #     if (!dw.error) dw.data.rep <- do.call(multiStageDynamics, arg.dw)
  #     if (!dw.error) result[sim, "dw"]  <- val.fn(dw.data.rep$summary$cumulative.event.time)
  #     result[sim, "time.dw"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
  #     arg.dw$policy <- NULL; gc()
  #     rm(dw.data.rep); gc()
  #   }
  #
  cat("\n******************************\n")
  cat ("6. Estimation - zero-order model for Simulation",sim, ":",generate_failure_method,"\n")
  if (!skip.zom) {
    if ("package:dtrSurv" %in% search()) {
      detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
    }
    cat ("  6. zero-order model - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    set.seed(train_seed + 5)
    optimal.zom <- do.call(itrSurv::itrSurv, c(arg.czmk2,
                                               list(mTry = 1,
                                                    nodeSize = 1e+9,
                                                    minEvent = 1)))
    # optimal.zom <- try(itrSurv:::itrSurv(data = data.df,
    #                                      txName = "Trt",
    #                                      models = models_itr,
    #                                      usePrevTime = TRUE,
    #                                      tau = tau,
    #                                      timePoints = "uni",
    #                                      nTimes = 100,
    #                                      criticalValue1 = criterion_phase1,
    #                                      criticalValue2 = criterion_phase2,
    #                                      evalTime = crit.value_phase1,
    #                                      splitRule = NULL,#ifelse(criterion == "mean", "mean", "logrank"),
    #                                      ERT = TRUE, uniformSplit = TRUE, replace = FALSE,
    #                                      randomSplit = 0.2,
    #                                      minEvent = 1, nodeSize = 1e+4, # zero-order model !!!!
    #                                      nTree = 300, mTry = 1,
    #                                      pooled = FALSE,
    #                                      tol1 = tol1,#c(0.1,0),
    #                                      stratifiedSplit = 0))
    zom.error <- class(optimal.zom)[1] == "try-error"
    arg.zom$policy <- if (!zom.error) optimal.zom
    policy_zom <- arg.zom$policy
    # print(arg.zom$policy@phaseResults[["EndPointPhase2Results"]]@optimal@optimalTx)
    if (!zom.error) result[sim, "zom_n_phase2"] <- mean(arg.zom$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) # result["zom_n_phase2"]
    rm(optimal.zom); gc()
    
    cat ("  6. zero-order model - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    if (!zom.error){
      arg_zom <<-arg.zom
      set.seed(train_seed + 10)
      zom.data.rep <- do.call(gdata_CR, arg.zom)
      rep_zom <<- zom.data.rep
      result[sim,"zom_survival"] = overall_survival_val.fn(zom.data.rep) #result["zom_survival"]
      result[sim,"zom_endpoint"] = endpoint_val.fn(zom.data.rep) #result["zom_endpoint"]
      # result[sim, "zom_percent.censor"] <- mean(zom.data.rep$status==0, na.rm = TRUE)
      # result[sim, paste0("zom_cause.", 1:n.causes)] <-
      #   sapply(1:n.causes, function(s) mean(zom.data.rep$status == s))
      
      # # within zom method, mean survival time from those who have action -1 and mean survival time for those who have action 1
      # A_mean = zom.data.rep %>% filter(action == -1) %>% overall_survival_val.fn()
      # B_mean = zom.data.rep %>% filter(action == 1) %>% overall_survival_val.fn()
      # A_mean_cr = zom.data.rep %>% filter(action == -1) %>% endpoint_val.fn()
      # B_mean_cr = zom.data.rep %>% filter(action == 1) %>% endpoint_val.fn()
      # trt_result[sim_count+4, c("method", "surv_A", "surv_B", "endpoint_A", "endpoint_B")] = c("zom", A_mean, B_mean, A_mean_cr, B_mean_cr)
      # rm(A_mean); rm(B_mean); rm(A_mean_cr); rm(B_mean_cr)
    }
    # if (!zom.error) result[sim, "zom_phase1"] <- val.fn_phase1(zom.data.rep$event.time)
    # if (!zom.error) result[sim, "zom_phase2"] <- val.fn_phase2(zom.data.rep$event.time)
    result[sim,"time.zom"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result[sim, "time.zom"]
    arg.zom$policy <- NULL; gc()
    rm(zom.data.rep); gc()
  }
  
  # comparing only those who continue to Phase 2 based on stopping rule:
  # for those who continue to P2, look at their CIF
  p2 = which(predd_surv_czmk_eval[["Stopping_Ind"]] == 0)
  if (!skip.czmk){
    P2_czmk <<- rep_czmk[p2, ]
    m1 = mean(P2_czmk$CIF_eval)
    m11 = mean(P2_czmk$OS_eval)
  } else{
    m1 = NA
    m11 = NA
  }
  if (!skip.csk){
    P2_csk <<- rep_csk[p2, ]
    m2 = mean(P2_csk$CIF_eval)
    m21 = mean(P2_csk$OS_eval)
  } else{
    m2 = NA
    m21 = NA
  }
  if (!skip.pmcr){
    P2_pmcr <<- rep_pmcr[p2, ]
    m3 = mean(P2_pmcr$CIF_eval)
    m31 = mean(P2_pmcr$OS_eval)
  } else{
    m3 = NA
    m31 = NA
  }
  if (!skip.aipwe){
    P2_aipwe <<- rep_aipwe[p2, ]
    m4 = mean(P2_aipwe$CIF_eval)
    m41 = mean(P2_aipwe$OS_eval)
  } else{
    m4 = NA
    m41 = NA
  }
  if (!skip.zom){
    P2_zom <<- rep_zom[p2, ]
    mz = mean(P2_zom$CIF_eval)
    mz1 = mean(P2_zom$OS_eval)
  } else{
    mz = NA
    mz1 = NA
  }
  P2_obs <<- rep_obs[p2, ]
  mo = mean(P2_obs$CIF_eval)
  mo1 = mean(P2_obs$OS_eval)
  
  p2.df = data.frame(methods = rep(NA, length(all_methods)))
  p2.df$methods = all_methods
  p2.df$OSeval = c(m11,m21,m31,m41,mz1,mo1)
  p2.df$CIFeval = c(m1,m2,m3,m4,mz,mo)
  
  
  # if (setting_seed == 1){
  #   arg.czmk$seed1 <- arg.zom$seed1 <- arg.obs$seed1 <- arg.obs.no.censor$seed1 <- NULL
  # }
  
  
  # realistic_train_truth <<- merge(train0.data, train1.data, by = "subj.id") %>%
  #   mutate(status = status.x,
  #          event.time = ifelse(event.time.x > event.time.y, event.time.x, event.time.y),
  #          action = ifelse(event.time.x > event.time.y, action.x, action.y)) %>%
  #   dplyr::select(-c(at.risk.x, at.risk.y, Z1.x,Z1.y,Z2.x,Z2.y,
  #                    # failure_t1.x, failure_t2.x, failure_t1.y, failure_t2.y,
  #                    status.x, status.y))
  # realistic_train_truth_list[[sim]] = realistic_train_truth
  #
  # train_truth <<- merge(train0.data, train1.data, by = "subj.id") %>%
  #   dplyr::select(-c(at.risk.x, at.risk.y, Z1.x,Z1.y,#Z2.x,Z2.y,
  #                    # failure_t1.x, failure_t2.x, failure_t1.y, failure_t2.y
  #                    )) %>%
  #   mutate(best.action = ifelse(event.time.x > event.time.y, action.x, action.y),
  #          best.status = ifelse(event.time.x > event.time.y, status.x, status.y),
  #          best.time = ifelse(event.time.x > event.time.y, event.time.x, event.time.y))  %>%
  #   dplyr::select(-c(action.x, action.y)) %>%
  #   mutate(obs.time0 = round(event.time.x,2),
  #          obs.time1 = round(event.time.y,2),
  #          status0 = status.x,
  #          status1 = status.y
  #   ) %>%
  #   dplyr::select(subj.id, best.time, best.action, best.status,
  #                 obs.time0, obs.time1, status0, status1)
  #
  # test_truth <<- merge(test0.data.rep, test1.data.rep, by = "subj.id") %>%
  #   dplyr::select(-c(Z1.x,Z1.y,Z2.x,Z2.y,
  #                    # failure_t1.x, failure_t2.x, failure_t1.y, failure_t2.y
  #                    )) %>%
  #   mutate(best.action.OS = ifelse(OS_eval.x > OS_eval.y, action.x, action.y),
  #          best.time.OS = ifelse(OS_eval.x > OS_eval.y, OS_eval.x, OS_eval.y),
  #          best.action.CIF = ifelse(CIF_eval.x < CIF_eval.y, action.x, action.y),
  #          best.time.CIF = ifelse(CIF_eval.x < CIF_eval.y, CIF_eval.x, CIF_eval.y))  %>%
  #   dplyr::select(-c(action.x, action.y)) %>%
  #   mutate(OS_eval0 = round(OS_eval.x,2),
  #          OS_eval1 = round(OS_eval.y,2),
  #          CIF_eval0 = round(CIF_eval.x,2),
  #          CIF_eval1 = round(CIF_eval.y,2),
  #          best.action = ifelse(best.action.OS == best.action.CIF, as.numeric(best.action.OS), 9999)
  #   ) %>%
  #   dplyr::select(subj.id, best.time.OS, best.action.OS, best.time.CIF, best.action.CIF,
  #                 OS_eval0, OS_eval1,
  #                 CIF_eval0, CIF_eval1,
  #                 best.action)
  #best.status, obs.time0, obs.time1, status0, status1)
  
  # View(cbind(true = test_truth$best.action, czmk = rep_czmk$action, csk = rep_csk$action, pmcr = rep_pmcr$action, zom = rep_zom$action, obs = rep_obs$action))
  #
  # realistic_test_truth <<- merge(test0.data.rep, test1.data.rep, by = "subj.id") %>%
  #   mutate(status = status.x,
  #          event.time = ifelse(event.time.x > event.time.y, event.time.x, event.time.y),
  #          action = ifelse(event.time.x > event.time.y, action.x, action.y)) %>%
  #   dplyr::select(-c(Z1.x,Z1.y,Z2.x,Z2.y,
  #                    # failure_t1.x, failure_t2.x, failure_t1.y, failure_t2.y,
  #                    status.x, status.y))
  # realistic_test_truth_list[[sim]] = realistic_test_truth
  
  
  # #OS
  # r00[sim, "sim"] = sim
  # r00[sim, "obs_training_os_trtprop"] = mean(train_truth$best.action == obs_1$action)
  # if (!skip.czmk) {r00[sim, "czmk_training_os_trtprop"] = mean(train_truth$best.action == policy_czmk@phaseResults[["FinalOptimalTx_Recc"]])}
  # if (!skip.csk) {r00[sim, "csk_training_os_trtprop"] = mean(train_truth$best.action == policy_csk@stageResults[[1]]@optimal@optimalTx)}
  # # if (!skip.pmcr) {r00[sim, "pmcr_training_os_trtprop"] = mean(train_truth$best.action == policy_csk@stageResults[[1]]@optimal@optimalTx)}
  # if (!skip.zom) {r00[sim, "zom_training_os_trtprop"] = mean(train_truth$best.action == policy_zom@phaseResults[["FinalOptimalTx_Recc"]])}
  # #CAUSE1
  # train_cause1_index = data.df %>% filter(status == 1) %>% dplyr::select(subj.id) %>% unlist() %>% as.vector()
  # train_truth1 = train_truth %>% filter(subj.id %in% train_cause1_index) %>% dplyr::select(-c(best.status, status0, status1, obs.time0, obs.time1))
  # r00[sim, "obs_training_cause1_trtprop"] = mean(train_truth1$best.action == obs_1$action[train_cause1_index])
  # if (!skip.czmk) {r00[sim, "czmk_training_cause1_trtprop"] = mean(train_truth1$best.action == policy_czmk@phaseResults[["FinalOptimalTx_Recc"]][train_cause1_index])}
  # if (!skip.csk) {r00[sim, "csk_training_cause1_trtprop"] = mean(train_truth1$best.action == policy_csk@stageResults[[1]]@optimal@optimalTx[train_cause1_index])}
  # if (!skip.zom) {r00[sim, "zom_training_cause1_trtprop"] = mean(train_truth1$best.action == policy_zom@phaseResults[["FinalOptimalTx_Recc"]][train_cause1_index])}
  #
  # #OS
  # r00[sim, "obs_testing_os_trtprop"] = mean(test_truth$best.action == rep_obs$action)
  # if (!skip.czmk) {r00[sim, "czmk_testing_os_trtprop"] = mean(test_truth$best.action == rep_czmk$action)}
  # if (!skip.csk) {r00[sim, "csk_testing_os_trtprop"] = mean(test_truth$best.action == rep_csk$action)}
  # if (!skip.zom) {r00[sim, "zom_testing_os_trtprop"] = mean(test_truth$best.action == rep_zom$action)}
  
  
  # print("sourcing SM01.")
  # source("SM01.True_ContinueP2_Ratios_Eval.R")
  # View(true2)
  
  ### saving and cleaning
  # print(90)
  # View(apply(result, 2, mean, na.rm = TRUE))
  # print(91)
  if (savingrds == TRUE){
    message('saving rds tmp')
    saveRDS(result, filename.tmp) # saving the temporary results
  }
  gc()
  cat("--- End of Simulation", sim, "---\n")
} # for (sim in 1:n.sim) but removed for parallelizing
# below is parallelization - uncomment if getting rid of for loop and sim index
# return(result)
# }

# below is parallelization
# # Number of cores to use for parallelization
# num_cores <- 5
# results_list = list()
# # Run the simulations in parallel
# results_list <- mclapply(1:n.sim, run_simulation, mc.cores = num_cores)


message("\n END OF SIMS \n")
# result
end_time = Sys.time()

message("End of Script: CR00_Simulation_Body.R.")
sprintf("Overall Time Took: %s", round(end_time - start_time,2))



# # Reshape the data into long format
# df_long <- melt(r00, id = "sim")
# # Extract method, A, and B from variable column
# df_long2 <- transform(df_long,
#                       method = sub("^(.*?)_.*$", "\\1", variable),
#                       data_type = sub("^.*?_(.*?)_.*$", "\\1", variable),
#                       endpoint_type = sub("^.*?_.*?_(.*?)_.*$", "\\1", variable),
#                       ignore = sub("^.*_([^_]+)$", "\\1", variable),
#                       proptrt = value)
#
# # Create the final data frame
# df_long3 <- df_long2[, c("sim","method", "data_type", "endpoint_type", "proptrt")]
# df_long3$method = factor(df_long3$method, levels = c("czmk", "csk", "zom", "obs"))
# df_long3$data_type = factor(df_long3$data_type, levels = c("training","testing"))
# df_long3$endpoint_type = factor(df_long3$endpoint_type, levels = c("os", "cause1"))
# # View(df_long3 %>% arrange(sim, endpoint_type, data_type,method))

