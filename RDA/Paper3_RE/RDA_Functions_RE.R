create_model_formulas <- function(outcome, covariates, formula = TRUE) {
  formula_string <- paste(paste(outcome), " ~", paste(covariates, collapse = " + "))
  if (formula == TRUE){
    formula_object <- as.formula(formula_string)
  } else{
    formula_object = formula_string
  }
  return(formula_object)
}

# single stage
weights_leuk <- function(data, lvls1) {
  # stage 1 was randomized
  propensity1 = sapply(lvls1, function(x) mean(data$A == x))
  propensity1 = sapply(data$A, function(s) propensity1[lvls1 == s])

  # IPCW
  Sc.hat1 <- coxph(Surv(T.1, 1 - d.1) ~ age, data = data) # Only 1 censored, so regress on only one var.
  Sc.hat1 <- exp( - predict(Sc.hat1, type = "expected"))
  Sc.hat <- Sc.hat1
  weight.censor = data$delta/Sc.hat
  return(list(propensity = propensity, weight.censor = weight.censor))
}

# used in weights_aric
St2 <- function(surv,
                time,
                t.eval,
                tau = Inf, exponential.tail = TRUE) {
  # tau is a placeholder for compatibility
  time = c(0, time)
  surv =  c(1, surv)
  max.t = max(time)
  min.s = min(surv)
  # print("--")
  # print(t.eval)
  # print(max.t)
  # print(min.s)
  if (t.eval > max.t) {
    if (exponential.tail) {
      S.t = 1 - (1-min.s)^(t.eval / max.t)
    } else {
      S.t = min.s
    }
  } else {
    # t.eval < max.t
    S.t = approxfun(time, surv)(t.eval)
  }
  return(S.t)
}

# data = testing_dataset; weight.formula.list = form.weight; covariates = covariate_names; event_indicator_string = event_indicator_string
weights_rda <- function(data, weight.formula.list, covariates, event_indicator_string) {
  require(dplyr)
  require(randomForest)
  require(randomForestSRC)

  # print("hi1")
  n = dim(data)[1]
  Tx.nm = sapply(weight.formula.list, function(x) x[2] %>%
                   as.character %>% gsub("(factor\\()(.*)(\\))", "\\2", .) )

  # print("hi2")
  # print(colnames(data))
  # print(weight.formula.list)
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
  # print("hi6")

  # IPCW
  # used to be Surv(time, 1 - d.1) --> make sure what we changed to is right (1-D.0)!
  # D.0 because overall survival and 1-D.0 because censoring survival probability
  rsf_censoring_string <- paste("Surv(TStop, 1 - ", event_indicator_string, ")", sep = "")
  # rfsrc requires all datatypes to NOT be character, so converting to factor
  data = data %>%
    mutate_if(is.character, as.factor)
  Sc.hat1 <<- rfsrc(create_model_formulas(rsf_censoring_string,
                                          covariates,
                                          formula = TRUE),
                    data = data)
  # print("hi7")
  # print(100)
  Sc.hat2 <<- sapply(1:dim(data)[1],
                     function(i)
                       St2(Sc.hat1$survival[i, ], #rows = people; cols = time of censored --> censoring survival probability for each person at each event time
                           Sc.hat1$time.interest, # the times that had someone censored in test set
                           t.eval = data$TStop[i])) # observed time for person i (evaluated time)
  # print("hi8")
  # clipping
  Sc.hat3 <<- Sc.hat2 %>% pmax(0.05) %>% pmin(0.95)
  weight.censor <<- data[[event_indicator_string]]/Sc.hat3 # overall survival
  # print("hi9")
  return(list(propensity = propensity2, weight.censor = weight.censor))
}

# test = arg.val$test
# actual = arg.val$actual
# estimated = get(opt_method)
# propensity = arg.val$propensity
# weight.censor = arg.val$weight.censor
# criterion = arg.val$criterion
# tau = arg.val$tau
getValue <- function(test_dat,
                     actual,
                     estimated,
                     propensity,
                     weight.censor,
                     criterion,
                     tau,
                     endpoint1) {
  # test_dat = testing_dataset_surv;
  # actual = test.tmp;
  # propensity = weight_prop %>% pull(weight); #weight$propensity[test_indices];
  # weight.censor = weight_censor %>% pull(weight); #weight$weight.censor[test_indices];
  # criterion = criterion;
  # tau = tau;
  # endpoint1 = endpoint1;
  # method = "ZOM"
  # opt_method <- paste0("opt.rule.", method)
  # estimated = get(opt_method)
  if (!(endpoint1 %in% c("survival", "mff"))) {
    stop("error: value must be either 'survival' or 'mff' for now")
  }
  
  # weight = apply(actual == estimated, 1, all, na.rm = TRUE) / propensity
  # print(1)
  all2 = function(x) {
    # if everything is NA, return NA. Otherwise, TRUE only if all is TRUE.
    # The naive all() returns NA even if there are only TRUEs except NAs.
    na.index = is.na(x)
    if (all(na.index)) return(NA)
    all(x[!na.index])
  }
  # print(2)
  weight_p =
    apply(actual == estimated, 1, all2) %>%
    # When there is at least one NA,
    # 1. all() does not return TRUE                     all(c(NA, NA, T)) = NA; all(c(T, T, T)) = TRUE
    # 2. all() returns NA if there is no FALSE          all(c(NA, NA, F, T)) = FALSE;  all(c(NA, NA)) = NA
    # 3. all() returns FALSE if there is at least one FALSE
    # When the second stage is not available in test set, the NA-match cases should still be counted. NA => 1.
    {ifelse(is.na(.), 1, as.numeric(.))} %>%
    "/" (propensity)
  if (endpoint1 == "survival"){
    weight.censor = weight.censor
    testY = test_dat[, "TStop"] #[, "obs_time"]
  } else{
    weight.censor = rep(1,length(weight.censor))
    testY = test_dat[, "recur"]/test_dat[, "TStop"]
  }
  weight = weight_p * weight.censor
  # evaluate(test_dat[, "T.0"], weight = weight, criterion = criterion, tau = tau)
  # print(3)
  if (criterion[1] == "mean" | criterion[1] == "area") {
    # print(3.5)
    mean(pmin(tau, testY) * weight)/ mean(weight)
  } else {
    mean(as.numeric(testY >= as.numeric(criterion[2])) * weight)/ mean(weight)
  }
  # print(4)
}

# ET <- function(surv, time, tau) {
#   tau.index = max(which(time <= tau))
#   time.last = time[tau.index]
#   surv = surv[1:tau.index]
#   time = time[1:tau.index]
#   surv = c(1, surv)
#   time.diff = c(time, tau) - c(0, time)
#   trunc.mean = sum(surv * time.diff)
#   trunc.mean
# }
# St <- function(surv, time, t.eval, tau) {
#   # tau is a placeholder for compatibility
#   t.eval = as.numeric(t.eval)
#   t.index.inf = max(which(time <= t.eval))
#   t.index.sup = min(which(time >= t.eval))
#   t.inf = time[t.index.inf]
#   t.sup = time[t.index.sup]
#   s.inf = surv[t.index.inf]
#   s.sup = surv[t.index.sup]
#   if (t.sup == t.inf) {
#     S.t = s.inf
#   } else {
#     S.t = s.inf -  (s.inf - s.sup) * (t.eval - t.inf)/(t.sup - t.inf)
#   }
#   S.t
# }
