test = test;
actual = test.tmp;
estimated = opt.rule.CSK;
# propensity = weight$propensity;
# weight.censor = weight$weight.censor;
criterion = value.criterion;
tau = tau

propensity = 1;
weight.censor = 0.1;

all2 = function(x) {
  # if everything is NA, return NA. Otherwise, TRUE only if all is TRUE.
  # The naive all() returns NA even if there are only TRUEs except NAs.
  na.index = is.na(x)
  if (all(na.index)) return(NA)
  all(x[!na.index])
}

weight = 
  apply(actual == estimated, 1, all2) %>%   
  # When there is at least one NA,
  # 1. all() does not return TRUE                     all(c(NA, NA, T)) = NA; all(c(T, T, T)) = TRUE
  # 2. all() returns NA if there is no FALSE          all(c(NA, NA, F, T)) = FALSE;  all(c(NA, NA)) = NA
  # 3. all() returns FALSE if there is at least one FALSE
  # When the second stage is not available in test set, the NA-match cases should still be counted. NA => 1.
  {ifelse(is.na(.), 1, as.numeric(.))} %>% 
  "/" (propensity)
weight = weight * weight.censor

testY = test[,"obs_time"]; #test[,"T.0"]

if (criterion[1] == "mean") {
  mean(pmin(tau, testY) * weight)/ mean(weight)
} else {
  mean(as.numeric(testY >= as.numeric(criterion[2])) * weight)/ mean(weight)
}
