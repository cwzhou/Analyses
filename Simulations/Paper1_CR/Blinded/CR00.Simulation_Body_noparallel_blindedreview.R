#### This script is run automatically from CR02.Simulation_Run.R. DO not run alone.

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
  mean(data$OS_eval, na.rm = TRUE)
}
endpoint_val.fn <- function(data) {
  # mean truncated years lost due to cause 1 (priorty cause) or
  mean(data$CIF_eval, na.rm = TRUE)
}
trt1.fn <- function(data) {
  # mean truncated years lost due to cause 1 (priorty cause) or
  mean(data$action==1, na.rm = TRUE)
}

if (criterion_phase1 %in% c("mean", "area", "surv.area")) {
  crit.eval.os = "mean"
  splitRule = "mean_surv"# for itrSurv Phase 1
} else if (criterion_phase1 %in% c("prob", "mean.prob.combo")) {
  crit.eval.os = "prob"
  splitRule = "logrank_surv" # for itrSurv Phase 1
}

# NOTE: below is assuming endpoint = CR aka we want mean LESS THAN critvalue
if (criterion_phase2 %in% c("mean", "area", "surv.area")) {
  crit.eval.cif = "mean"
} else if (criterion_phase2 %in% c("prob", "mean.prob.combo")) {
  crit.eval.cif = "prob"
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
      paste(method, "trt1", sep = "_"),
      paste("time", method, sep = ".")
    )
  } else {
    c(
      paste(method, "survival", sep = "_"),
      paste(method, "endpoint", sep = "_"),
      paste(method, "trt1", sep = "_"),
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
  
  if (sensitivity == 1){
    ztype1 = 0 #mix of uniform and normal
  } else{
    ztype1 = 2 #normal distribution
  }
  
  arg_list =   list(
    N = n, tau = tau, # structural parameters
    ztype = ztype1, 
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

### simulation
r00 = data.frame(sim = n.sim) #this is for sim for loop: data.frame(sim = n.sim)
r00[column_names1] <- NA

# Create the data frame with the 'sim' column
result <- data.frame(sim = n.sim_start:n.sim_end) # this is for parallelizaitn: data.frame(sim.no = sim) #this is for sim for loop: data.frame(sim = 1:n.sim)
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

true2 = true3 = list()
print(Sys.time())
for (sim in n.sim_start:n.sim_end){
  
  row_index = sim - n.sim_start + 1
  
  cat("\n\n#################################")
  cat("\n######### Simulation ",sim, "#########")
  cat("\nEndpoint: ",endpoint)
  if (endpoint == "CR"){
    cat("\nNumber of CR: ", M)
  }
  cat("\n#################################\n")

  train_seed = sim*10000 + init_seed*3
  test_seed = train_seed + 30306
  
  cat ("1. Data for Simulation",sim, ":",generate_failure_method,"\n")
  tt(1)
  
  # obs.data
  set.seed(init_seed*sim + 505)
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
  
  message("data.df")
  data.df <<- obs.data %>%
    mutate(D.0 = ifelse(status> 0,1,0), # event from any cause indicator
           D.1 = ifelse(status==1,1,0), # event from cause 1 indicator
           D.2 = ifelse(status==2,1,0), # event from cause 2 indicator
           obs_time = event.time,
           Trt = action) %>%
    dplyr::select(obs_time, Trt, status, D.0, D.1, D.2, contains("Z")) %>%
    arrange(obs_time) # MUST SORT ASCENDING ORDER FOR OBS_TIME
  
  # create unique observed failure times for survival phase 1
  timePointsSurvival = data.df %>%
    filter(D.0 == 1) %>%
    dplyr::select(obs_time) %>%
    unlist(use.names = FALSE) %>%
    unique()
  # for CR only
  timePointsEndpoint = timePointsSurvival # we can do this bc CR is subset of OS failure times
  # create unique observed endpoint times for endpoint phase 2
  
  if (!skip.pmcr){
    data.df.pmcr <<- data.df %>%
      dplyr::select(obs_time, Trt, status, starts_with("Z")) %>%
      mutate(Trt = ifelse(Trt == -1, 0, Trt)) # Trt has to be 0/1
  } else{
    message("skipping pmcr data manipulation because skip.pmcr = TRUE")
  }
  
  if (!skip.aipwe){
    data.df.aipwe <<- aipwe_data_format(data.df)
  } else{
    message("skipping aipwe data manipulation because skip.aipwe = TRUE")
  }
  
  result[row_index, "training_percent.censor"] <- mean(data.df$status==0, na.rm = TRUE) #result["training_percent.censor"]
  # result[paste0("training_cause.", 1:n.causes)]
  result[row_index, paste0("training_cause.", 1:n.causes)] <-
    sapply(1:n.causes, function(s) mean(data.df$status == s))
  
  # # # obs policy value
  message("using test_seed to generate obs testing data")
  set.seed(test_seed)
  obs.data.rep <- do.call(gdata_CR, arg.obs.no.censor) # no censoring for eval sets
  rep_obs <<- obs.data.rep
  
  result[row_index, "obs_survival"] <- overall_survival_val.fn(obs.data.rep)#val.fn_phase1(obs.data.rep$event.time)
  result[row_index, "obs_endpoint"] <- endpoint_val.fn(obs.data.rep)
  result[row_index, "obs_trt1"] <- trt1.fn(obs.data.rep)
  
  result[row_index, "time.obs"] <- tt(2, reset = TRUE, unit = "min")["elapsed"] #result["time.obs"]
  # # print(flowchart(obs.data$output)) # if i fix this in 02_Simulation_Functions.R for CR then I can use, otherwise delete.
  # # transforming data from an array format to a data.frame format
  # # data.df = output2observable(obs.data)
  rm(obs.data.rep); gc()
  
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
                                                     minEventEnd = mindeath, # minimum number of subjects with events
                                                     minEventSurv = mindeath, # minimum number of subjects with events
                                                     nodeSizeEnd = nodesize,
                                                     nodeSizeSurv = nodesize # this is needed only for endpoint=RE.
                                                )))
    print(table(optimal.czmk@phaseResults[["FinalOptimalTx_Recc"]]))
    czmk.error <- class(optimal.czmk)[1] == "try-error"
    arg.czmk$policy <- if (!czmk.error) optimal.czmk
    policy_czmk <<- arg.czmk$policy
    #training
    predd_surv_czmk <<- predd_surv
    predd_ep_czmk <<- predd_ep
    # print(arg.czmk$policy@phaseResults[["EndPointPhase2Results"]]@optimal@optimalTx)
    if (!czmk.error) result[row_index, "czmk_n_phase2"] <- mean(arg.czmk$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) #result["czmk_n_phase2"]
    rm(optimal.czmk); gc()
    
    cat ("  \n 2. czmk - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_czmk <<- arg.czmk
    if (!czmk.error) {
      set.seed(train_seed + 2)
      czmk.data.rep <- do.call(gdata_CR, arg.czmk)
      predd_surv_czmk_eval <<- predd_surv
      predd_ep_czmk_eval <<- predd_ep
      rep_czmk <<- czmk.data.rep
      
      indiv_St_opt <<- as.data.frame(docalling1[["optimal"]]@optimalY[,])
      indiv_St_trt0 <<- docalling1[["value"]][["Func"]][[1]]
      indiv_St_trt1 <<- docalling1[["value"]][["Func"]][[2]]
      # message("checking priority-cause CIF curve for first 3 timepoints for everyone")
      indiv_CIF_opt <<- as.data.frame(docalling2[["optimal"]]@optimalY[,])
      indiv_CIF_trt0 <<- docalling2[["value"]][["Func"]][[1]]
      indiv_CIF_trt1 <<- docalling2[["value"]][["Func"]][[2]]
      # plotting_indiv_curves_after_prediction.R

      result[row_index, "czmk_survival"] = overall_survival_val.fn(czmk.data.rep) #result["czmk_survival"]
      result[row_index, "czmk_endpoint"] = endpoint_val.fn(czmk.data.rep) #result["czmk_endpoint"]
      result[row_index,"czmk_trt1"] = trt1.fn(czmk.data.rep)
    }
    
    result[row_index, "time.czmk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.czmk"]
    arg.czmk$policy <- NULL; gc()
    rm(czmk.data.rep); gc()
  }
  
  cat("\n******************************\n")
  # estimation
  cat("3. csk - Cho et al for Simulation",sim, ":",generate_failure_method,"\n")
  n.stages = 1
  if (!skip.csk) {
    if ("package:itrSurv" %in% search()) {
      detach("package:itrSurv", unload = TRUE, character.only = TRUE)
    }
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
    set.seed(train_seed + 3)
    optimal.csk <- do.call(dtrSurv::dtrSurv, c(arg.csk2, list(nodeSize = nodesize, minEvent = mindeath )))
    csk.error <- class(optimal.csk)[1] == "try-error"
    arg.csk$policy <- if (!csk.error) optimal.csk
    policy_csk <<- arg.csk$policy
    rm(optimal.csk); gc()
    
    cat ("  \n 3. csk - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_csk <<- arg.csk
    
    if (!csk.error) {
      set.seed(test_seed)
      csk.data.rep <- do.call(gdata_CR, arg.csk)
      rep_csk <<- csk.data.rep
      result[row_index, "csk_survival"] = overall_survival_val.fn(csk.data.rep) #result["csk_survival"]
      result[row_index, "csk_endpoint"] = endpoint_val.fn(csk.data.rep) # result["csk_endpoint"]
      result[row_index,"csk_trt1"] = trt1.fn(csk.data.rep)
     }
    
    result[row_index, "time.csk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.csk"]
    arg.csk$policy <- NULL; gc()
    rm(csk.data.rep); gc()
    if ("package:dtrSurv" %in% search()) {
      detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
    }
    library(itrSurv)
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
    set.seed(train_seed + 4)
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
    set.seed(train_seed + 4)
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
      set.seed(test_seed)
      pmcr.data.rep <- do.call(gdata_CR, arg.pmcr)
      rep_pmcr <<- pmcr.data.rep
      result[row_index, "pmcr_survival"] = overall_survival_val.fn(pmcr.data.rep) #result["pmcr_survival"]
      result[row_index, "pmcr_endpoint"] = endpoint_val.fn(pmcr.data.rep) #result["pmcr_endpoint"]
      result[row_index,"pmcr_trt1"] = trt1.fn(pmcr.data.rep)
      }
    
    result[row_index, "time.pmcr"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.pmcr"]
    arg.pmcr$policy <- NULL; gc()
    rm(pmcr.data.rep); gc()
    }
  
  
  cat("\n******************************\n")
  # estimation
  cat ("5. aipwe for Simulation",sim, ":",generate_failure_method,"\n")
  if (!skip.aipwe) {
    cat ("  5. aipwe - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    # models_aipwe <- paste0("Trt ~ ",
    # paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
    arg.aipwe2 = list(data = data.df.aipwe, # b/c trt is 0/1
                      pp.v = ncov,#c(2,1,2),#c(1,1,1,2), #ncov,
                      tau1 = t0_aipwe, # t0 needs to be tuned
                      tune = c(0.001,0.01,0.5,1,
                               seq(0.1,400,
                                   length.out=16))
    )
    set.seed(train_seed + 5)
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
    
    cat ("  \n 5. aipwe - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    arg_aipwe <<- arg.aipwe
    if (!aipwe.error) {
      set.seed(test_seed)
      aipwe.data.rep <- do.call(gdata_CR, arg.aipwe)
      rep_aipwe <<- aipwe.data.rep
      result[row_index, "aipwe_survival"] = overall_survival_val.fn(aipwe.data.rep) #result["aipwe_survival"]
      result[row_index, "aipwe_endpoint"] = endpoint_val.fn(aipwe.data.rep) #result["aipwe_endpoint"]
      result[row_index,"aipwe_trt1"] = trt1.fn(aipwe.data.rep)
    }
    
    result[row_index, "time.aipwe"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result["time.aipwe"]
    arg.aipwe$policy <- NULL; gc()
    rm(aipwe.data.rep); gc()
  }
  
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
                                                    minEventEnd = 1, # minimum number of subjects with events
                                                    minEventSurv = 1, # minimum number of subjects with events
                                                    nodeSizeEnd = 1e9,
                                                    nodeSizeSurv = 1e9 # this is needed only for endpoint=RE.
                                               )))
    print(table(optimal.zom@phaseResults[["FinalOptimalTx_Recc"]]))

    zom.error <- class(optimal.zom)[1] == "try-error"
    arg.zom$policy <- if (!zom.error) optimal.zom
    policy_zom <- arg.zom$policy
    # print(arg.zom$policy@phaseResults[["EndPointPhase2Results"]]@optimal@optimalTx)
    if (!zom.error) result[row_index, "zom_n_phase2"] <- mean(arg.zom$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) # result["zom_n_phase2"]
    rm(optimal.zom); gc()
    
    cat ("  6. zero-order model - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
    if (!zom.error){
      arg_zom <<-arg.zom
      set.seed(test_seed)
      zom.data.rep <- do.call(gdata_CR, arg.zom)
      rep_zom <<- zom.data.rep
      result[row_index,"zom_survival"] = overall_survival_val.fn(zom.data.rep) #result["zom_survival"]
      result[row_index,"zom_endpoint"] = endpoint_val.fn(zom.data.rep) #result["zom_endpoint"]
      result[row_index,"zom_trt1"] = trt1.fn(zom.data.rep)
    }
    result[row_index,"time.zom"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result[row_index, "time.zom"]
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
  
  ### saving and cleaning
  if (savingrds == TRUE){
    message('saving rds tmp')
    saveRDS(result, filename.tmp) # saving the temporary results
  }
  message("The result for up to simulation ", sim)
  print(result)
  gc()
  cat("--- End of Simulation", sim, "---\n")
} # for (sim in 1:n.sim) 

message("\n END OF SIMS \n")
# result
end_time = Sys.time()

message("End of Script: CR00_Simulation_Body.R.")
sprintf("Overall Time Took: %s", round(end_time - start_time,2))

