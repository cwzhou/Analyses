set.seed(2023)
N = 1000
u_10 = runif(N)
for (alpha in 4*c(0.0625,0.125, 0.25)){
  print(sprintf("#### alpha: %s #### ", alpha))
  tt = gumbelbiexp(u_10, alpha, 12.1)
}

s1 <- survfit(Surv(tt, rep(1,N)) ~ 1)
s = s1$surv
LAM1 = -log(s)[1:N-1]
time1 = s1$time[1:N-1]
df = data.frame(time1,LAM1)
ggplot(data = df, aes(x = time1, y = LAM1)) +
  geom_line() +
  xlab("Time") +  # Customize x-axis label
  ylab("LAM1") +  # Customize y-axis label
  ggtitle("Line Plot of LAM1 over Time: tt")  # Add a title to the plot

samples = rexp(N,1)
s1 <- survfit(Surv(samples, rep(1,N)) ~ 1)
s = s1$surv
LAM1 = -log(s)[1:N-1]
time1 = s1$time[1:N-1]
df = data.frame(time1,LAM1)
ggplot(data = df, aes(x = time1, y = LAM1)) +
  geom_line() +
  xlab("Time") +  # Customize x-axis label
  ylab("LAM1") +  # Customize y-axis label
  ggtitle("Line Plot of LAM1 over Time: rexp")  # Add a title to the plot



s1 <- survfit(Surv(samples, rep(1,N)) ~ 1)
df1 = data.frame(LAM1 = s1$surv, time1 = s1$time)
ggplot(data = df, aes(x = time1, y = LAM1)) +
  geom_line() +
  xlab("Time") +  # Customize x-axis label
  ylab("LAM1") +  # Customize y-axis label
  ggtitle("Line Plot of LAM1 over Time")  # Add a title to the plot

