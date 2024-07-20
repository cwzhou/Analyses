train.tmp = train %>% transmute(Trt = train[, Tx.nm])

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
  training_dataset = train %>%
    dplyr::select("obs_time",
                  "Trt",
                  any_of(covariate_names),
                  event_indicator_string)

  # weights
  weight_name <- paste0("weight.", suffix1)
  weight = weights_rda(data = training_dataset,
                       weight.formula.list = form.weight,
                       covariates = covariate_names,
                       event_indicator_string = event_indicator_string)
  arg.val_name <- paste0("arg.val.", suffix1)
  arg.val = list(test = training_dataset,
                 actual = train.tmp,
                 propensity = weight$propensity,
                 weight.censor = weight$weight.censor,
                 criterion = criterion,
                 tau = tau)

  for (method in loop_methods) {
    err_method <- paste0("err.", method)
    opt_method <- paste0("opt.rule.", method)
    if (!get(err_method)) {

      if (method == "CZMK"){
        estimated1 = CZMK.i@phaseResults[["FinalOptimalTx_Recc"]]
      } else if (method == "ZOM"){
        estimated1 = ZOM.i@phaseResults[["FinalOptimalTx_Recc"]]
      } else if (method == "CSK"){
        estimated1 = CSK.i@stageResults[[1]]@optimal@optimalTx
      } else if (method == "CSKzom"){
        estimated1 = CSKzom.i@stageResults[[1]]@optimal@optimalTx
      }


      # Truncated Mean/Prob (anything other than OS/PC is counted as censored)
      train_eval_result_list_tmp[[method]][suffix1] = do.call(getValue,
                                                              c(arg.val,
                                                                list(estimated = estimated1)))
    }
  }

  # observed
  obs1 <- paste0("observed.", suffix1)
  train_eval_result_list_tmp[["observed"]][suffix1] = do.call(getValue,
                                                              c(arg.val[-which(names(arg.val) == "propensity")],
                                                                list(estimated = train.tmp,
                                                                     propensity = 1)))
}

train_eval_result_list[[cv]] = train_eval_result_list_tmp
