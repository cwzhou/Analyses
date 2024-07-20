
df_cov = data.frame(c(1:N), z, A)
colnames(df_cov) = c("ID", colnames(z), "Trt")

# generating censoring time
if (ctype == 0){cc <- rexp(N,cparam)}
if (ctype == 1){
  if(is.null(censor_max)){censor_max = tau}
  if(is.null(censor_min)){censor_min = 0}
  cc <- runif(N,min=censor_min,max=censor_max)}
# print(sprintf("censoring time cc: %s",cc))

# generating failure time from cause 1 based on Fine-Gray paper (1999):
u1 <- runif(N)  # Observed value
failure_t1 = failure_t2 = obs_time_failureCR = numeric(0)
for (i in 1:N){
  u1i = u1[i]
  pred.hazard1i = pred.hazard1[i]
  # Use backsolve_t function
  failure_t1[i] <- backsolve_t1(u1i, mass_p, pred.hazard1i)
}
failure_t2_rate = exp(pred.hazard2)
failure_t2 = rexp(N, rate = failure_t2_rate)
pred.hazard2 <<- pred.hazard2
failre_t2_rate <<- failure_t2_rate
failre_t2 <<- failure_t2


df_times0 = data.frame(ID = 1:N,
                       Time_Failure1 = failure_t1,
                       Time_Failure2 = failure_t2,
                       Time_Censor = cc,
                       Time_Tau = tau
)
df_times0 <<-df_times0
# print(length(df_times0$Time_Failure2));print(length(df_times0$Time_Failure1));
# print(length(df_times0$Time_Censor));print(length(df_times0$Time_Tau))
df_times = df_times0 %>%
  mutate(obs_time_failureCR = pmin(Time_Failure1, Time_Failure2),
         obs_time = pmin(Time_Failure1,Time_Failure2,Time_Censor,Time_Tau),
         status = ifelse(obs_time == Time_Failure1, 1,
                         ifelse(obs_time == Time_Failure2, 2,
                                0))
  )

df_surv1 = df_times %>%
  # overall death indicator (failure from any cause = 1; 0 o/w)
  mutate(indD = ifelse(obs_time_failureCR <= pmin(Time_Censor, Time_Tau), 1, 0))
df_surv1$obs_time = df_surv1$obs_time*365
df_surv2 = inner_join(df_cov, df_surv1, by = "ID")
