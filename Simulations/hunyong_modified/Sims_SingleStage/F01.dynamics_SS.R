library(MASS); library(ggplot2); library(dplyr); library(Rcpp)

#' @return event.time observed time (min of T1,T2,C)
#' @return T1 competing failure time 1
#' @return T2 competing failure time 2
#' @return event_cause = 1 if event.time = T1; = 2 if event.time = T2, and = 0 if event.time = C
#' @return delta 1(T_tilde <= C); C = censoring time (overall survival indicator)

# time.max = 1; tick = tick;
# cause1_prob = cause1_prob;
# censor_rate = 0.2;
# pred.hazard = 0;
# full_data = ifelse(hidden_data, 2, 0)

gdata_CR <- function(time.max = 3,
                     tick = 0.01,
                     cause1_prob = 0.5,
                     censor_rate = 0.2,
                     pred.hazard = 0,
                     multiplier = 2,
                     full_data = FALSE) {
  # full_data: present all hidden event times (failure, treatment, censoring) and cause-specific failure
  #            if (full_data == 2) present all hidden event times without cause-specific failure
  # outputs :
  # event.time = X = min(T1,T2,C), delta = 1(not censored), event_cause = 1(failure cause) and 0 = censored

  timeframe = seq(0, time.max, by = tick)
  n.time = length(timeframe)

  ## event times
  # Following the simulation setup from Fine-Gray (1999) paper
  # Subdistribution for type1 failures given by unit exponential mixture with mass 1-caues1_prob at inf when Zi = 0 and uses proportional
  # subdistribution for nonzero covariates:
  # Pr(Ti \le t, type = 1|Zi) = 1-[1-cause1_prob{1-exp(-t)}]^exp(predHazardFn(action, covariate, cause = 1))
  # Subdistribution for type2 failures obtained by taking Pr(type=2|Zi)=1-Pr(type=1|Zi) and using exponential
  # distribution with rate exp(predHazardFn(action, covariate, cause = 2) for Pr(Ti \le t|type=2,Zi)
  # Censoring times generated from uniform[a,b] distribution

  # Cause1 failure time
  u21 = runif(1)
  hazard = exp(pred.hazard)



  # censoring
  C = runif(N, min = a, max = b)

  data = data.frame(ID = 1:N,
                    failure1_time = failure_t1,
                    failure2_time = failure_t2,
                    obs_failureCR_time = pmin(failure_t1, failure_t2),
                    censoring_time = C
  )

  # Calculate obs_time as the minimum of failure times and censoring time
  data$obs_time <- do.call(pmin, c(data[c("failure1_time",
                                          "failure2_time",
                                          "censoring_time")], na.rm = TRUE))

  # Create status variable based on obs_time and different failure times
  data$status <- ifelse(data$obs_time == data$failure1_time, 1,
                        ifelse(data$obs_time == data$failure2_time, 2,
                               0))

  type1FailureRate = sum(data$status == 1)/N
  type2FailureRate = sum(data$status == 2)/N
  censoringRate = sum(data$status == 0)/N
  message("type 1 failure rate: ", type1FailureRate)
  message("type 2 failure rate: ", type2FailureRate)
  message("censoring rate: ", censoringRate)
  View(data)
  data = data %>% dplyr::select(ID, obs_time, status)






  surv = exp(-cumsum(hazard))
  failure.time = which(surv <= runif(1))[1] * tick # NA if administratively censored
  if (is.na(failure.time)) failure.time <- Inf

  trt.time = Inf

  ## summary statistics
  X = min(failure.time, trt.time)
  gamma = failure.time <= trt.time
  x.index = which(timeframe == X)
  # censoring
  censor.time = rexp(1, censor_rate)
  XX = min(X, censor.time)
  delta = (X <= censor.time)
  print(delta)

  if (!delta) { #if censored (when delta = 0)
    # print("!delta")
    # message("delta: ",delta)
    output = c(event.time = XX, gamma = NA, delta = delta)
    # message("full_data:",full_data)
    if (full_data == 1) {
      return(list(statistics = output,
                  times = c(failure.time = failure.time, treatment.time = trt.time,
                            censor.time = censor.time),
                  surv = surv))
    } else if (full_data == 2) {
      return(c(output, failure.time = failure.time, treatment.time = trt.time,
               censor.time = censor.time))
    } else {
      return(output)
    }
  }


  # cat(failure.time, " ", trt.time, " ",  gamma,"\n")
  output = c(event.time = X, gamma = gamma, delta = delta)

  if (full_data == 1) {
    return(list(statistics = output,
                times = c(failure.time = failure.time, treatment.time = trt.time,
                          censor.time = censor.time),
                surv = surv))
  } else if (full_data == 2) {
    return(c(output, failure.time = failure.time, treatment.time = trt.time,
             censor.time = censor.time))
  } else {
    return(output)
  }
}
gdata_CR.vec <- Vectorize(gdata_CR, vectorize.args = c("pred.hazard"))

#cum_ou in Rcpp
cppFunction(
  "NumericVector cum_ou(NumericVector vec, double reverting_const) {
    const int n = vec.size();
    NumericVector y(n);
    y[0] = vec[0];
    for (int i=1; i < n; ++i) {
    y[i] = y[i] + y[i-1] * (1 - reverting_const);
    }
    return y;
    }")
