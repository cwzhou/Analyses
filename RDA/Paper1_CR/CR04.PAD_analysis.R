# NOTE: CURRENTLY AIPWE TRAIN/TEST DATASETS DON'T USE MULTI-LEVEL FACTORS!!! 

library(caret);
library(purrr);
library(dplyr);
library(survival)
library(DTRreg);
library(randomForestSRC);
#library(dtrSurv) # we only use library for CSK method if-statements
library(itrSurv);library(ggplot2); library(tidyverse);
library(MASS); library(dplyr);
local = 1
setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR/")
source("F01.Simulation_Functions.R") # includes F02.ComparatorMethod_Functions.R
setwd("../../RDA/Paper1_CR/")
source("F03.RDA_Functions.R")    # All the library and source files
############################################################
######### temporary testing things out parameters ##########
############################################################
criterion_phase1 = "mean"
rule2 = "gray_cr" # other option: csh_cr
############################################################
################# methods for comparison ###################
############################################################
# NOTE: CURRENTLY AIPWE TRAIN/TEST DATASETS DON'T USE MULTI-LEVEL FACTORS!!! NEED TO FIX!!!
rda_methods = c("CZMK", "CSK", "PMCR", "AIPWE",
                "ZOM", "CSKzom", "observed")
skip_method <- !c(TRUE,TRUE,TRUE,TRUE,
                 TRUE,!TRUE,TRUE);
assign_skip_function(rda_methods, skip_method)
skipped_methods <- rda_methods[skip_method]
loop_methods <- rda_methods[!rda_methods %in% c(skipped_methods, "observed")]

############################################################
#################### 0. parameters ####################
############################################################
tol1_param = c(0.1,0,0.3,0.01)
t0_crit = 365/2 #2200 # 6-mo survival
pooled1 = FALSE # stratified = lower nodesize; pooled = can have larger nodesize
tau = 365
K = 10 #300 # number of CV
endpoint = "CR" # endpoint
Tx.nm = "Trt"
timepoints = seq(0, sqrt(tau), length.out = 1000)^2

## other parameters
priority_cause = 1 # death
nodesize = 5
mindeath = round(sqrt(c(nodesize)), 0)
Ntree = 300
ert = TRUE;
rs = 0.2 # randomSplit = 0.2

############################################################
###################### UNC PAD Cohort ######################
############################################################
pad_df = read.csv("/Volumes/McGinigle_Kate_IRB18-1153/Analysis/pad_df.csv",
                  stringsAsFactors = TRUE)
dataset_name = "pad"
# code for reading in data
data = pad_df %>%
  # make sure Treatment is "Trt" and time is "obs_time"
  mutate(Trt = ifelse(Trt == "Endovascular", 1, 0)) %>%
  # replace any covariate names with space or - with _
  rename_with(~ gsub(" ", "_", .x), contains(" ")) %>%
  rename_with(~ gsub("-", "_", .x), contains("-"))

# get survialtimes
s.times = data %>%
  filter(D.0 == 1) %>%
  dplyr::select(obs_time) %>%
  arrange(obs_time) %>%
  unlist() %>% as.numeric()
# get CR priority cause times
pc.times = data %>%
  filter(D.1 == 1) %>%
  dplyr::select(obs_time) %>%
  arrange(obs_time) %>%
  unlist() %>% as.numeric()

# Get all column names of the dataframe
all_column_names <- colnames(data)
# Define the column names you want to exclude
excluded_column_names <- c("ID", "Trt", "obs_time",
                           "status", "D.0", "D.1", "D.2")
# Exclude the specified column names
covariate_names <- setdiff(all_column_names, excluded_column_names)
# covar_names = dat0 %>% dplyr::select(-c(D.0, D.1, D.2, status, Trt, obs_time)) %>% names()
# paste0(paste(covar_names, collapse = ", "))

# manipulating dataset to obtain truncated status,D.0,D.1,D.2
dat0 = data %>%
  # mutate(status_new = ifelse(obs_time > tau, 0, status)) %>%
  # mutate(D.0_new = ifelse(status > 0, 1, 0),
  #        D.1_new = ifelse(status == 1, 1, 0),
  #        D.2_new = ifelse(status == 2, 1, 0)) %>%
  dplyr::select(obs_time, Trt,
                status, D.0, D.1, D.2,
                # status_new, D.0_new, D.1_new, D.2_new,
                !!!covariate_names)
# head(dat0);range(dat0$obs_time)

  #################### 1. data preprocessing
  # sorting dataset by obs_time (ascending response time)
  dat = as.data.frame(dat0) %>% arrange(obs_time) # required for CR
  # dat = dat0 %>%
  #   mutate(time = pmin(time, tau),
  #      last.fu = ifelse(time == tau, 1, last.fu),
  #      d.1 = ifelse(time == tau & last.fu == 1, 1, d.1),
  #      d.2 = ifelse(time == tau, ifelse(last.fu == 1, NA, 1), d.2),
  #      delta = ifelse(time == tau, 1, delta))

  mean(dat$status==0, na.rm = TRUE)

  lvls = lapply(Tx.nm, function(x) levels(dat[, x]))
  for (x in 1:length(lvls)) {if (is.null(lvls[[x]])) lvls[[x]] = 0:1}


  ############################################################
  ###################### 0.2 criterion #######################
  ############################################################
  # for (criterion_phase1 in c("area", "mean", "mean.prob.combo")){ #for (criterion_phase1 in c("mean", "mean.prob.combo")){ # for (dtr_criterion in c("mean", "surv.mean")) {
  print(criterion_phase1)
  criterion_phase2 = criterion_phase1
  if (criterion_phase1 == "mean.prob.combo"){
    dtr_criterion = "surv.mean"
  } else if (criterion_phase1 == "prob"){
    dtr_criterion = "surv.prob"
  } else{
    dtr_criterion = "mean"
  }
  if (!(criterion_phase1[1] %in% c("mean", "area"))) {
    criterion_phase1[2] = t0_crit
    dtr_criterion[2] = criterion_phase1[2]
    rule1 = "logrank_surv"
  } else {
    rule1 = "mean_surv"
    criterion_phase1[2] = NA
  }
  # rule1 =  ifelse(criterion_phase1[1] == "mean" | criterion_phase1[1] == "area", "mean_surv", "logrank_surv")

  if (!(criterion_phase2[1] %in% c("mean", "area"))) {
    criterion_phase2[2] = t0_crit
  } else {
    criterion_phase2[2] = NA
  }
  if (is.null(criterion_phase1[2]) | is.na(criterion_phase1[2])){
    t0_pmcr = tau/2 # change this !!
    t0_aipwe = tau/2
    crit_tmp = ""
  } else{
    t0_pmcr = criterion_phase1[2]
    t0_aipwe = criterion_phase1[2]
    crit_tmp = round(as.numeric(criterion_phase1[2]), 1)
  }

  cat("\ncriterion_phase1 =",criterion_phase1[1],
      "\nvalue.os =", criterion_phase1[2],
      "\nvalue.cif =", criterion_phase2[2],
      "\ntau =", tau, "\n")

  #################### Overall Survival PM models
  modelOS = create_model_formulas("Surv(obs_time, D.0)",
                                  covariate_names,
                                  formula = FALSE)
  form.OS <- modelOS %>% as.formula()

  #################### propensity models
  modelPr_PMCR = create_model_formulas("Trt", #trt must be 0/1
                                       covariate_names,
                                       formula = TRUE)
  modelPr = create_model_formulas("factor(Trt)",
                                  covariate_names,
                                  formula = FALSE)
  form.weight <- list(modelPr %>% as.formula)

  ## cross-validation
  K = K
  set.seed(2024, kind="Mersenne-Twister")
  cv.insample  = createDataPartition(dat[, Tx.nm],
                                     p = 0.8, # aipwe requires a certain number of people to work (aka can't train on anything too little sample size)
                                     list = FALSE,
                                     times = K)
  cv.outsample = map_dfc(1:K, ~which(!(1:dim(dat)[1] %in% cv.insample[, .]))) %>% as.matrix

  # skeleton
  values_colsnames = c(paste0(rda_methods, ".OS"),
                       paste0(rda_methods, ".PC"),
                       "ns.CZMK", "ns.CSK",
                       "train_cens", "train_cause1", "train_cause2")#
  values <- matrix(NA, K, length(values_colsnames),
                   dimnames = list(1:K, values_colsnames))
  train_eval_result_list = train_eval_result_list_tmp = list()

  nm = paste0(criterion_phase1[1],
              crit_tmp, "_rule1",
              rule1, "_rule2", rule2, "_tau", tau)
  fnm = paste0("./3_output/", endpoint, "/",dataset_name,"/", Sys.Date(), "/") # folder name
  if (!dir.exists(sprintf("./3_output/%s", endpoint))) dir.create(sprintf("./3_output/%s", endpoint))
  if (!dir.exists(sprintf("./3_output/%s/%s", endpoint, dataset_name))) dir.create(sprintf("./3_output/%s/%s", endpoint, dataset_name))
  if (!dir.exists(sprintf("./3_figure/%s", endpoint))) dir.create(sprintf("./3_figure/%s", endpoint))
  if (!dir.exists(sprintf("./3_figure/%s/%s", endpoint, dataset_name))) dir.create(sprintf("./3_figure/%s/%s", endpoint, dataset_name))
  if (!dir.exists(fnm)) dir.create(fnm)
  rds = paste0("_", Sys.Date(), ".rds")

  tab = data.frame(s1 = rep(NA, K),
                   # s2 = NA,
                   total = NA)

  for (cv in 1:K) {
    cat(cv, "th cv.\n")
    set.seed(cv)
    in.cv = cv.insample[, cv]   # insample index
    # print(tail(in.cv))
    out.cv = cv.outsample[, cv] # outsample index
    train = dat[in.cv, ] %>% as.data.frame()
    test  = dat[out.cv,] %>% as.data.frame()
    if (skip.PMCR != TRUE){
      train_pmcr = train %>% dplyr::select(-c(D.0, D.1, D.2)) %>%
        mutate(across(everything(), convert_factors_to_numeric))
    }
    if (skip.AIPWE != TRUE){
      train_aipwe0 = train %>%
        dplyr::select(-c("ischemia", "woundClass", "maxRutherfordClass", "Race"))
      train_aipwe = aipwe_data_format(train_aipwe0)
      test_aipwe = test %>%
        dplyr::select(-c("ischemia", "woundClass", "maxRutherfordClass", "Race"))
      # account for factors
    }
    values[cv, "train_cens"] = mean(train$status==0); # proportion of people censored
    values[cv, "train_cause1"] = mean(train$status==1, na.rm = T)
    values[cv, "train_cause2"] = mean(train$status==2, na.rm = T)

    ############################################################################################################
    #### A. rule estimation ####################################################################################
    ############################################################################################################

    ### A1. The proposed method
    if (skip.CZMK != TRUE){
      print(sprintf("Running CZMK for %s CV", cv))
      priority_vector <- c(0, priority_cause)

      # Use lapply to create a list of sprintf statements
      models_itr <- lapply(priority_vector, function(x) {
        paste0(sprintf("Surv(obs_time, D.%s) ~ ",
                       ifelse(x == 0, "0", as.character(x))),
               gsub(".*~\\s*", "", modelOS)) %>% as.formula
        })
      args.CZMK <- list(data = train,
                        endPoint = "CR",
                        epName = "status",
                        yName = "obs_time",
                        txName = Tx.nm,
                        models = models_itr,
                        tau = tau,
                        # timePoints = timepoints,
                        timePointsSurvival = s.times,
                        timePointsEndpoint = s.times,
                        criticalValue1 = criterion_phase1[1],
                        criticalValue2 = criterion_phase2[1],
                        evalTime = as.numeric(criterion_phase1[2]),
                        splitRule1 = rule1,
                        splitRule2 = rule2,
                        ERT = ert, uniformSplit = ert, replace = !ert,
                        randomSplit = rs, nTree = Ntree, mTry = 6,
                        pooled = pooled1,
                        tol1 = if (criterion_phase1[1] == 'mean.prob.combo') tol1_param
                        else tol1_param[1:2],#c(0.1,0), #c(0.1,0),
                        stratifiedSplit = FALSE)# actual fitting
      values[cv, "ns.CZMK"] = nodesize
      set.seed(cv)
      CZMK.i <-
        try(do.call(itrSurv::itrSurv,
                    c(args.CZMK, list(nodeSizeSurv = nodesize,
                                      nodeSizeEnd = nodesize,
                                      minEventSurv = mindeath,
                                      minEventEnd = mindeath))))
      err.CZMK = class(CZMK.i)[1] == "try-error"
      } else{
        err.CZMK = TRUE
        }

    ### A2. Cho et al (2022)
    if (skip.CSK != TRUE){
      library(dtrSurv)
      print(sprintf("Running CSK for %s CV", cv))
      args.CSK <- list(data = train,
                       txName = Tx.nm,
                       models = form.OS,
                       usePrevTime = FALSE, tau = tau, timePoints = timepoints,
                       criticalValue = dtr_criterion[1],
                       evalTime = as.numeric(dtr_criterion[2]),
                       splitRule = ifelse(dtr_criterion[1] == "mean", "mean", "logrank"),
                       ERT = ert, uniformSplit = ert, replace = !ert,
                       randomSplit = rs, nTree = Ntree, mTry = 6,
                       pooled = pooled1, stratifiedSplit = FALSE)
      # actual fitting
      values[cv, "ns.CSK"] = nodesize
      set.seed(cv)
      CSK.i <- try(do.call(dtrSurv::dtrSurv,
                           c(args.CSK, list(nodeSize = nodesize,
                                            minEvent = mindeath))))
      err.CSK = class(CSK.i)[1] == "try-error"
      if ("package:dtrSurv" %in% search()) {
        detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
      }
      } else{
        err.CSK = TRUE
        }

    ### A3. PMCR - 2021
    if (skip.PMCR != TRUE){
      print(sprintf("Running PMCR for %s CV", cv))
      message("unrestricted regime")
      args.PMCR <- args.PMCR2 <-
        list(Time="obs_time",
             Event="status",
             formula=modelPr_PMCR, #trt needs to be 0/1
             data=train_pmcr,
             rgenoud=FALSE,
             Restrict=FALSE,
             propscore="logistic",
             t0=t0_pmcr)
      set.seed(cv)
      temp.unrestr.fit<-try(do.call(PMCR, args.PMCR))
      err.PMCRunrestr = class(temp.unrestr.fit)[1] == "try-error"
      if (!err.PMCRunrestr){
        #range of alpha # alpha needs to be tuned
        alps<-c(temp.unrestr.fit$Fbeta2[2],
                seq(round(temp.unrestr.fit$Fbeta2[2],2)+0.01,
                    round(temp.unrestr.fit$Fbeta1[2],2)+0.03,
                    0.01))
        alp<-alps[2] # alp needs to be tuned based on alps
        message("restricted regime")
        M_pmcr<-1e+5
        args.PMCR2$Restrict = TRUE
        args.PMCR2$M = M_pmcr
        args.PMCR2$alp = alp
        temp.restr.fit<-try(do.call(PMCR, args.PMCR2))
        # temp.restr.fit<-PMCR(Time="obs_time",
        #                      Event="status", # CR event indicator
        #                      formula=modelPr_PMCR, # Trt must be 0/1
        #                      data=train_pmcr,
        #                      rgenoud=FALSE,
        #                      Restrict=TRUE,
        #                      propscore="logistic",
        #                      t0=t0_pmcr,
        #                      alp=alp,
        #                      M=M_pmcr)
        # dat1 = temp.restr.fit$data
        # cov = dat1 %>%
        #   dplyr::select(-c(id, obs_time, status, pix.hat)) %>%
        #   mutate(int = rep(1, nrow(dat1))) %>%
        #   dplyr::select(int, covariate_names)
        PMCR.i <- temp.restr.fit$beta3
        err.PMCR = class(PMCR.i)[1] == "try-error"
        } else{
          err.PMCR = TRUE
          }
      } else{
        err.PMCR = TRUE
        }

    ### A4. AIPWE
    if (skip.AIPWE != TRUE){

      ### A4. AIPWE - 2022
      print(sprintf("Running AIPWE for %s CV", cv))

      args.AIPWE = list(data_list = train_aipwe,
                        pp.v = (ncol(train_aipwe$Z1)-1)/2, #minus trt; divided by 2 b/c of interactions
                        tau1 = as.numeric(t0_aipwe),
                        tune = c(0.001,0.01,0.5,1,
                                 seq(0.1,400,
                                     length.out=16)))
      AIPWE.i <-
        try(do.call(aipwe.fit,
                    c(args.AIPWE)))
      err.AIPWE = class(AIPWE.i)[1] == "try-error"

      # aif = aipwe.fit(data_list = train_aipwe,
      #                 pp.v = (ncol(train_aipwe$Z1)-1)/2, #minus trt; divided by 2 b/c of interactions
      #                 tau1 = as.numeric(t0_aipwe),
      #                 tune = c(0.001,0.01,0.5,1,
      #                          seq(0.1,400,
      #                              length.out=16))) #eta0-eta_{ncov}
      # AIPWE.i <- aif

    } else{
      err.AIPWE = TRUE
    }

    ### A6. zero-order model
    if (skip.ZOM != TRUE){
      print(sprintf("Running ZOM for %s CV", cv))
    ZOM.i <-
      try(do.call(itrSurv::itrSurv, c(args.CZMK,
                             list(nodeSize = 1e+4,
                                  minEvent = 1e+4))))
    err.ZOM = class(ZOM.i)[1] == "try-error"

    if (!err.ZOM) {
      zom.pred =
        data.frame(
          ZOM.i@phaseResults$FinalOptimalTx_Recc[1])
      tab[cv, 1] = ZOM.i@phaseResults$FinalOptimalTx_Recc[1]
      cat("ZOM freq table\n"); tab[, 1] %>% table %>% print
    }
    } else{
      err.ZOM = TRUE
    }


    ### A7. CSKzero-order model
    if (skip.CSKzom != TRUE){
      library(dtrSurv)
      print(sprintf("Running CSKzom for %s CV", cv))
    CSKzom.i <-
      try(do.call(dtrSurv::dtrSurv, c(args.CSK,
                             list(nodeSize = 1e+4,
                                  minEvent = 1e+4))))
    err.CSKzom = class(CSKzom.i)[1] == "try-error"

    if (!err.CSKzom) {
      CSKzom.pred =
        data.frame(
          CSKzom.i@stageResults[[1]]@optimal@optimalTx[1])
      tab[cv, 1] = CSKzom.i@stageResults[[1]]@optimal@optimalTx[1]
      cat("CSKzom freq table\n");
      tab[, 1] %>% table %>% print
    }
    if ("package:dtrSurv" %in% search()) {
      detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
    }
    } else{
      err.CSKzom = TRUE
    }

    message("Applying rules to test set.")
    ############################################################################################################
    #### B. Rules applied to a test set ########################################################################
    ############################################################################################################
    ### B0. skeletons for the predicted Trt's
    opt.rule.pred = data.frame(Trt = rep(NA, dim(cv.outsample)[1]))
    opt.rule.pred[,1] = factor(NA, levels = lvls[[1]])

    for(method in loop_methods) {
      assign(paste("opt.rule.", method, sep = ""),
             opt.rule.pred,
             envir = .GlobalEnv)
    }
    # opt.rule.CZMK <- opt.rule.CSK <- opt.rule.PMCR <- opt.rule.ZOM <- opt.rule.CSKzom <- opt.rule.pred
    # opt.rule.GKRF <- opt.rule.GKLM <-
    elig = !test[, Tx.nm] %>% is.na
    ### B1. the proposed method
    if (!err.CZMK) {
      for (phase in 1:2){
        opt.CZMK_phase =
          predict(CZMK.i,
                  newdata = test[elig, ],
                  Phase = phase,
                  endPoint = "CR")
        if (phase == 1){
          action1 <<- opt.CZMK_phase$optimal@optimalTx
          StopatP1 <<- opt.CZMK_phase$optimal@Ratio_Stopping_Ind
        }
        if (phase ==2){
          action2 <<- opt.CZMK_phase$optimal@optimalTx
        }
      }
      tmp_act = as.data.frame(cbind(P1=action1, P2=action2, StopatP1))
      tmp_act = tmp_act %>%
        mutate(Final = ifelse(StopatP1 == 1,
                              P1,
                              P2))
      opt.CZMK = tmp_act$Final
      opt.rule.CZMK[elig, ] = factor(opt.CZMK,
                                     levels = lvls[[1]]) # %>% as.numeric()
      print(table(opt.rule.CZMK, useNA = "always"))
    }
    ### B2. Cho et al (2022)
    if (!err.CSK) {
      library(dtrSurv)
      opt.CSK =
        predict(CSK.i,
                # newdata = test[elig, ],
                newdata = test[elig, ],
                stage = 1)
      opt.rule.CSK[elig, ] = factor(opt.CSK$optimal@optimalTx,
                                    levels = lvls[[1]]) # %>% as.numeric()
      print(table(opt.rule.CSK, useNA = "always"))
      if ("package:dtrSurv" %in% search()) {
        detach("package:dtrSurv", unload = TRUE, character.only = TRUE)
      }
    }

    ### B3. PMCR - 2021
    if (!err.PMCR) {

      terms1 <- terms(modelPr_PMCR, data = test[elig, ])
      cov2 <- model.matrix(terms1,test[elig, ]) %>%
        as.matrix()
      opt.PMCR = ifelse(
        (cov2 %*% as.matrix(PMCR.i)) > 0,
        1, 0) %>% as.vector()

      opt.rule.PMCR[elig, ] = factor(opt.PMCR,
                                     levels = lvls[[1]]) # %>% as.numeric()
      print(table(opt.rule.PMCR, useNA = "always"))
    }

    ### B4. AIPWE - 2022
    if (!err.AIPWE) {

      cov2 <- test_aipwe[elig,] %>%
        dplyr::select(-c(obs_time,Trt,status,D.0,D.1,D.2)) %>%
        as.matrix()
      eta0 = AIPWE.i[1] #
      eta = AIPWE.i[-1]
      opt.AIPWE =
        ifelse((eta0 + cov2 %*% as.matrix(eta)) < 0, 1, 0) %>%
        as.vector()
      opt.rule.AIPWE[elig, ] = factor(opt.AIPWE,
                                     levels = lvls[[1]]) # %>% as.numeric()
      print(table(opt.rule.AIPWE, useNA = "always"))
    }

    ### B6. zero-order model
    if (!err.ZOM) {
      opt.rule.ZOM[elig, ] = zom.pred[1, ]
      opt.rule.ZOM[elig, ] = factor(zom.pred[1, ], levels = lvls[[1]]) # %>% as.numeric()
      print(table(opt.rule.ZOM, useNA = "always"))
    }

    ### B7. CSKzero-order model
    if (!err.CSKzom) {
      opt.rule.CSKzom[elig, ] = CSKzom.pred[1, 1]
      opt.rule.CSKzom[elig, ] = factor(CSKzom.pred[1, 1],
                                       levels = lvls[[1]]) # %>% as.numeric()
      print(table(opt.rule.CSKzom, useNA = "always"))
    }

    ############################################################################################################
    #### C. Value estimation ###################################################################################
    ############################################################################################################
    test.tmp = test %>% transmute(Trt = test[, Tx.nm])

    if (criterion_phase1[1] != criterion_phase2[1]){
      message("WARNING: CRITERION_PHASE1 AND CRITERION_PHASE2 ARE DIFFERENT - MAKE SURE THAT IS WHAT YOU WANT")
    }
    suffixes = c("OS", "PC")
    for (suffix in 1:length(suffixes)){
      suffix1 = suffixes[suffix]
      if (suffix1 == "OS"){
        event_indicator_string = "D.0"
        criterion = criterion_phase1
      } else if (suffix1 == "PC"){
        event_indicator_string = toString(paste0("D.",priority_cause))
        criterion = criterion_phase2
      }
      testing_dataset = test %>%
        dplyr::select("obs_time",
                      "Trt",
                      any_of(covariate_names),
                      all_of(event_indicator_string))

      # weights
      weight_name <- paste0("weight.", suffix1)
      weight = weights_rda(data = testing_dataset,
                           weight.formula.list = form.weight,
                           covariates = covariate_names,
                           event_indicator_string = event_indicator_string)
      arg.val_name <- paste0("arg.val.", suffix1)
      arg.val = list(test = testing_dataset,
                        actual = test.tmp,
                        propensity = weight$propensity,
                        weight.censor = weight$weight.censor,
                        criterion = criterion,
                        tau = tau)

      for (method in loop_methods) {
        err_method <- paste0("err.", method)
        opt_method <- paste0("opt.rule.", method)
        if (!get(err_method)) {
          # Truncated Mean/Prob (anything other than OS/PC is counted as censored)
          values[cv, paste0(method, ".", suffix1)] <-
            do.call(getValue,
                    c(arg.val,
                      list(estimated =
                             get(opt_method))))
        }
      }

      # observed
      obs1 <- paste0("observed.", suffix1)
      values[cv, obs1] <- do.call(getValue,
                                  c(arg.val[-which(names(arg.val) == "propensity")],
                                    list(estimated = test.tmp,
                                         propensity = 1)))
    } # end of suffix

    source("Training_Eval_Results.R")

    print(values[cv, ])
    print(apply(values, 2, mean, na.rm = TRUE)) # across columns (over all cv iterations)

    ############################################################################################################
    #### C. saving the results #################################################################################
    ############################################################################################################
    attr(values, "spec") = data.frame(criterion1 = criterion_phase1[1],
                                      criterion2 = criterion_phase2[1],
                                      criterion1.s = as.numeric(criterion_phase1[2]),
                                      criterion2.s = as.numeric(criterion_phase2[2]),
                                      rule1 = rule1,
                                      rule2 = rule2,
                                      tau = tau,
                                      ert = ert,
                                      rs = rs,
                                      ntree = Ntree,
                                      nodesize = nodesize,
                                      mindeath = mindeath)
    # if (cv==1) saveRDS(CZMK.i, paste0(fnm, "czmk.cv",cv, "_", dataset_name, "_", nm, "_", cv, rds))
    # if (cv==1) saveRDS(CSK.i, paste0(fnm, "csk.cv", cv, "_", dataset_name, "_", nm, "_", cv, rds))
    # if (cv==1) saveRDS(PMCR.i, paste0(fnm, "pmcr.cv", cv, "_", dataset_name,"_", nm, "_", cv, rds))
    saveRDS(values, paste0(fnm, "values_cv",cv,"_values_", nm, rds))
  } # end of cv
  saveRDS(values, paste0(fnm,"Values_", K,"CV_", nm, rds))
  saveRDS(tab, paste0(fnm, "tab_ZOM_", nm, rds))
  # View(values)
# } # end of criterion.value_phase1
