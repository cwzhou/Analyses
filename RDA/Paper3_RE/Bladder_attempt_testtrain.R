# Assuming sub_ids and dat are already defined
K <- 5  # Number of cross-validation iterations
values <- data.frame(cv = 1:K, train_cens = numeric(K), test_cens = numeric(K))
cat(round(mean(dat_surv$Status_D == 1),3)*100, "% patients die.")
set.seed(2025)
for (cv in 1:K) {
  # Randomly sample 30% of IDs for the test set
  set.seed(cv*K)  # For reproducibility
  test_ids <- sample(sub_ids, size = floor(0.3 * length(sub_ids)))
  train_ids <- setdiff(sub_ids, test_ids)
  
  # Find the indices in 'dat' corresponding to test and train IDs
  test_indices <- which(dat$id %in% test_ids)
  train_indices <- which(dat$id %in% train_ids)
  
  # Subset the data for train and test sets
  train <- dat %>% 
    filter(id %in% train_ids) %>% 
    dplyr::select(-c(recur)) %>% 
    as.data.frame()
  
  test <- dat %>% 
    filter(id %in% test_ids) %>% 
    as.data.frame()
  
  # Debug print statements
  print(paste("CV Iteration:", cv))
  print(paste("Train size:", length(unique(train$id)), "Unique IDs with", nrow(train), "records"))
  print(paste("Test size:", length(unique(test$id)), "Unique IDs with", nrow(test), "records"))
  
  # Calculate the proportion of censored data conditionally on Status == 0
  train_censored <- train %>% filter(Status == 0)
  test_censored <- test %>% filter(Status == 0)
  
  values[cv, "train_cens"] <- mean(train_censored$Status_D == 0)
  values[cv, "test_cens"] <- mean(test_censored$Status_D == 0)
  cat(round(mean(train_censored$Status_D == 1),3)*100, "% patients in train set die.")
  cat(round(mean(test_censored$Status_D == 1),3)*100, "% patients in test set die.")
  
  model1 = "Surv(TStop, Status_D) ~ number + size" %>% as.formula()
  model2 = "Surv(TStart, TStop, Status) ~ number + size" %>% as.formula()
  models_RE = list(model1, model2)
  
  data_surv = train %>% filter(Status == 0)
  
  args.CZMK <- list(data = train,
                    endPoint = "RE",
                    idName = "id",
                    epName = "Status",
                    txName = "A",
                    models = models_RE,
                    timePointsSurvival = timePointsSurvival,
                    timePointsEndpoint = timePointsEndpoint,
                    tau = tau,
                    criticalValue1 = "mean",
                    criticalValue2 = "mean",
                    evalTime = 1,
                    splitRule1 = "mean_surv",
                    splitRule2 = "gray_re",
                    ERT = FALSE,
                    uniformSplit = TRUE,
                    replace = FALSE,
                    randomSplit = 0.2,
                    nTree = Ntree,
                    pooled = FALSE,
                    tol1 = tol1_param,
                    stratifiedSplit = 0.1)
  
  ############################################################################################################
  #### A. rule estimation ####################################################################################
  ############################################################################################################
  ### A1. The proposed method
  if (skip.CZMK != TRUE){
    print(sprintf("Running CZMK for %s CV", cv))
    
    values[cv, "ns1.CZMK"] = nodeSizeSurv
    values[cv, "ns2.CZMK"] = nodeSizeEnd
    set.seed(cv)
    CZMK.i <-
      try(do.call(itrSurv::itrSurv,
                  c(args.CZMK, list(nodeSizeSurv = nodeSizeSurv,
                                    nodeSizeEnd = nodeSizeEnd,
                                    minEventSurv = minEventSurv,
                                    minEventEnd = minEventEnd))))
    # values[cv, "CZMK.train.terminal"] = CZMK.i@value[["V1"]][["Et_survival"]]
    values[cv, "CZMK.PropPhase2"] = CZMK.i@value[["V2"]][["PropPhase2"]]
    # values[cv, "CZMK.train.RE"] = CZMK.i@value[["V3"]][["Et_mff"]]
    err.CZMK = class(CZMK.i)[1] == "try-error"
  } else{
    err.CZMK = TRUE
  }
  
  ### A6. zero-order model
  if (skip.ZOM != TRUE){
    print(sprintf("Running ZOM for %s CV", cv))
    set.seed(cv)
    ZOM.i <-
      try(do.call(itrSurv::itrSurv, c(args.CZMK,
                                      list(nodeSizeEnd = 1e+4,
                                           nodeSizeSurv = 1e+4,
                                           minEventEnd = 1e+4,
                                           minEventSurv = 1e+4))))
    # values[cv, "ZOM.train.terminal"] = ZOM.i@value[["V1"]][["Et_survival"]]
    values[cv, "ZOM.PropPhase2"] = ZOM.i@value[["V2"]][["PropPhase2"]]
    # values[cv, "ZOM.train.RE"] = ZOM.i@value[["V3"]][["Et_mff"]]
    err.ZOM = class(ZOM.i)[1] == "try-error"
    
    if (!err.ZOM) {
      zom.pred =
        data.frame(
          ZOM.i@phaseResults$FinalOptimalTx_Recc[1])
      tab[cv, 1] = ZOM.i@phaseResults$FinalOptimalTx_Recc[1]
      # cat("ZOM freq table\n"); tab[, 1] %>% table %>% print
    }
  } else{
    err.ZOM = TRUE
  }
  
  ############################################################################################################
  #### B. Rules applied to a test set ########################################################################
  ############################################################################################################
  message("Applying rules to test set.")
  ### B0. skeletons for the predicted Trt's
  opt.rule.pred = data.frame(Trt = rep(NA, length(unique(test$id))))
  opt.rule.pred[,1] = factor(NA, levels = lvls[[1]])
  
  for(method in loop_methods) {
    assign(paste("opt.rule.", method, sep = ""),
           opt.rule.pred,
           envir = .GlobalEnv)
  }
  # Create eligibility for each row based on non-missing values in `Tx.nm`
  elig0 <- test %>%
    dplyr::select(id, !!sym(Tx.nm)) %>%
    group_by(id) %>%
    dplyr::mutate(eligibility = all(!is.na(.data[[Tx.nm]]))) %>% # Check if all values are TRUE
    dplyr::mutate(eligibility = ifelse(any(!eligibility), FALSE, eligibility)) %>% # Set all to FALSE if any are FALSE
    ungroup()
  elig = elig0 %>%
    pull(eligibility)
  id_elig = elig0 %>%
    dplyr::select(id, eligibility) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    pull(eligibility)
  
  ### B1. the proposed method
  if (!err.CZMK) {
    for (phase in 1:2){
      opt.CZMK_phase =
        itrSurv::predict(CZMK.i,
                         newdata = test[elig, ],
                         Phase = phase,
                         epName1 = "Status",
                         endPoint = "RE")
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
    opt.rule.CZMK[id_elig, ] = factor(opt.CZMK,
                                      levels = lvls[[1]]) # %>% as.numeric()
    print(table(opt.rule.CZMK, useNA = "always"))
  }
  
  ### B6. zero-order model
  if (!err.ZOM) {
    opt.rule.ZOM[id_elig, ] = rep(zom.pred[1, ], length(id_elig))
    opt.rule.ZOM[id_elig, ] = factor(zom.pred[1, ], levels = lvls[[1]]) # %>% as.numeric()
    print(table(opt.rule.ZOM, useNA = "always"))
  }
  
  ############################################################################################################
  #### C. Value estimation ###################################################################################
  ############################################################################################################
  test.surv = test %>%
    filter(Status == 0)
  test.tmp = test.surv%>%
    transmute(Trt = test.surv[, Tx.nm])
  
  if (criterion_phase1[1] != criterion_phase2[1]){
    message("WARNING: CRITERION_PHASE1 AND CRITERION_PHASE2 ARE DIFFERENT - MAKE SURE THAT IS WHAT YOU WANT")
  }
  suffixes = c(".test.terminal", ".test.RE")
  for (suffix in 1:length(suffixes)){
    suffix1 = suffixes[suffix]
    if (suffix1 == ".test.terminal"){
      event_indicator_string = "Status_D"
      criterion = criterion_phase1
      endpoint1 = "survival"
    } else if (suffix1 == ".test.RE"){
      event_indicator_string = "Status"
      criterion = criterion_phase2
      endpoint1 = "mff"
    }
    testing_dataset = test %>%
      dplyr::select("TStart", "TStop",
                    "recur",
                    Tx.nm,
                    any_of(covariate_names),
                    all_of(event_indicator_string))
    testing_dataset_surv = test %>%
      filter(Status == 0) %>%
      dplyr::select("TStart", "TStop", "recur",
                    Tx.nm,
                    any_of(covariate_names),
                    all_of(event_indicator_string))
    
    # # weights
    # weight_name <- paste0("weight", suffix1)
    # weight = weights_rda(data = testing_dataset,
    #                      weight.formula.list = form.weight,
    #                      covariates = covariate_names,
    #                      event_indicator_string = event_indicator_string)
    # arg.val_name <- paste0("arg.val.", suffix1)
    # arg.val = list(test = testing_dataset,
    #                actual = test.tmp,
    #                propensity = weight$propensity,
    #                weight.censor = weight$weight.censor,
    #                criterion = criterion,
    #                tau = tau)
    
    weight_prop = cbind(id = test_ids,
                        weight = weight0$propensity[test_indices]) %>%
      as.data.frame()
    weight_censor = cbind(id = test_ids,
                          weight = weight0$weight.censor[test_indices]) %>%
      as.data.frame()
    # Expand weight to match the number of records for each person in the test dataset
    weight_prop_long <- test %>%
      # Join the weight dataset to the test dataset by 'id'
      left_join(weight_prop, by = "id") %>%
      pull(weight)
    weight_censor_long <- test %>%
      # Join the weight dataset to the test dataset by 'id'
      left_join(weight_censor, by = "id") %>%
      pull(weight)
    
    arg.val = list(test_dat = testing_dataset_surv,
                   actual = test.tmp,
                   propensity = weight_prop %>%
                     pull(weight), #weight$propensity[test_indices],
                   weight.censor = weight_censor %>%
                     pull(weight), #weight$weight.censor[test_indices],
                   criterion = criterion,
                   tau = tau,
                   endpoint1 = endpoint1)
    for (method in loop_methods) {
      err_method <- paste0("err.", method)
      opt_method <- paste0("opt.rule.", method)
      if (!get(err_method)) {
        # est1 = cbind(id = test_ids,
        #              est = get(opt_method)) %>% as.data.frame()
        # est2 = test %>%
        #   dplyr::select(id) %>%
        #   left_join(est1, by = "id") %>%
        #   pull(Trt) %>%
        #   as.data.frame()
        est2 = get(opt_method)
        print(head(est2))
        
        # Truncated Mean/Prob
        values[cv, paste0(method, suffix1)] <-
          do.call(getValue,
                  c(arg.val,
                    list(estimated = est2)))
        #get(opt_method))))
      }
    }
    # if (!err.CSK)  values[cv, "CSK"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.CSK)))
    # if (!err.zom)  values[cv, "ZOM"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.zom)))
    
    # observed
    obs1 <- paste0("OBS", suffix1)
    values[cv, obs1] <- do.call(getValue, c(arg.val, list(estimated = test.tmp)))
  } # end of suffix
  
  # source("Training_Eval_Results.R")
  message("Printing values")
  print(values[cv, ])
  message("Printing mean across cv iterations")
  print(apply(values, 2, mean, na.rm = TRUE)) # across columns (over all cv iterations)
  # View(values)
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
                                    ntree = Ntree)
  # saveRDS(values, paste0(fnm, "values_cv",cv,"_values_", nm, rds))
} # end of LOOCV

print(values)
cat("averaging across iterations")
print(apply(values, 2, mean, na.rm = TRUE)) # across columns (over all cv iterations)

