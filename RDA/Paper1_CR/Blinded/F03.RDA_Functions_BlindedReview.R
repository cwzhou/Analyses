create_model_formulas <- function(outcome, covariates, formula = TRUE) {
  formula_string <- paste(paste(outcome), " ~", paste(covariates, collapse = " + "))
  if (formula == TRUE){
    formula_object <- as.formula(formula_string)
  } else{
    formula_object = formula_string
  }
  return(formula_object)
}

weights_rda <- function(data, weight.formula.list, covariates, event_indicator_string) {
  require(dplyr)
  require(randomForest)
  require(randomForestSRC)

  # print("hi1")
  n = dim(data)[1]
  Tx.nm = sapply(weight.formula.list, function(x) x[2] %>%
                   as.character %>% gsub("(factor\\()(.*)(\\))", "\\2", .) )
  
  # propensity scores are estimated by RF using test set
  prop.model = randomForest(weight.formula.list[[1]],
                             data = data, )
  # print("hi4")
  propensity = predict(prop.model,
                        data, type = "prob")   # n x n.levels matrix
  # print("hi5")
  if (is.null(dim(propensity))) {
    # in case it is binary
    propensity = data.frame(1 - propensity, propensity)
    names(propensity) = 0:1
  }
  missing1 = is.na(propensity[, 1]) & !is.na(data[[event_indicator_string]]) # make sure D.0 is right. used to be d.1 (stage 1)             # Bookkeeping covariate missing subjects (neither censored nor failed)
  propensity1 = propensity[cbind(rownames(data),
                                  data[, Tx.nm[1]] %>%
                                    as.character)]  # A vector of predicted pi(A|X) with NA's
  propensity2 = ifelse(is.na(propensity1), 1, propensity1)             # Replace NA's with 1's.

  # IPCW
  # used to be Surv(time, 1 - d.1) --> make sure what we changed to is right (1-D.0)!
  # D.0 because overall survival and 1-D.0 because censoring survival probability
  rsf_censoring_string <- paste("Surv(obs_time, 1 - ", event_indicator_string, ")", sep = "")
  # rfsrc requires all datatypes to NOT be character, so converting to factor
  data = data %>%
    mutate_if(is.character, as.factor)
  Sc.hat1 <<- rfsrc(create_model_formulas(rsf_censoring_string,
                                           covariates,
                                           formula = TRUE),
                    data = data)
  Sc.hat2 <<- sapply(1:dim(data)[1],
                    function(i)
                      St2(Sc.hat1$survival[i, ], #rows = people; cols = time of censored --> censoring survival probability for each person at each event time
                          Sc.hat1$time.interest, # the times that had someone censored in test set
                          t.eval = data$obs_time[i])) # observed time for person i (evaluated time)
  # clipping
  Sc.hat3 <<- Sc.hat2 %>% pmax(0.05) %>% pmin(0.95)
  weight.censor <<- data[[event_indicator_string]]/Sc.hat3 # overall survival
  return(list(propensity = propensity2, weight.censor = weight.censor))
}

getValue <- function(test,
                     actual,
                     estimated,
                     propensity,
                     weight.censor,
                     criterion,
                     tau) {
  
  all2 = function(x) {
    # if everything is NA, return NA. Otherwise, TRUE only if all is TRUE.
    # The naive all() returns NA even if there are only TRUEs except NAs.
    na.index = is.na(x)
    if (all(na.index)) return(NA)
    all(x[!na.index])
  }

  weight_p =
    apply(actual == estimated, 1, all2) %>%
    # When there is at least one NA,
    # 1. all() does not return TRUE                     all(c(NA, NA, T)) = NA; all(c(T, T, T)) = TRUE
    # 2. all() returns NA if there is no FALSE          all(c(NA, NA, F, T)) = FALSE;  all(c(NA, NA)) = NA
    # 3. all() returns FALSE if there is at least one FALSE
    # When the second stage is not available in test set, the NA-match cases should still be counted. NA => 1.
    {ifelse(is.na(.), 1, as.numeric(.))} %>%
    "/" (propensity)
  weight = weight_p * weight.censor
  # evaluate(test[, "T.0"], weight = weight, criterion = criterion, tau = tau)
  testY = test[, "obs_time"]
  if (criterion[1] == "mean" | criterion[1] == "area") {
    mean(pmin(tau, testY) * weight)/ mean(weight)
  } else {
    mean(as.numeric(testY >= as.numeric(criterion[2])) * weight)/ mean(weight)
  }
}
