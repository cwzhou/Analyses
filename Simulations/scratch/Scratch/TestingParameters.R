
# function for generating failure time from Gumbel bivariate exponential distribution (1960)
gumbelbiexp <- function(u, alpha, y) {
  cc = 1-2*exp(-y); 
  print("1-2exp(-y)")
  # print(cc)
  rootsqd = (u-1)/(cc*alpha) + ((1+cc*alpha)/(2*cc*alpha))^2; #print(rootsqd)
  root = sqrt(rootsqd);
  # print("root")
  # print(root)
  expnegx_plus = (1+cc*alpha)/(2*cc*alpha) + root;
  expnegx_minus = (1+cc*alpha)/(2*cc*alpha) - root;
  print((expnegx_plus > 1))
  print(expnegx_minus)
  print((expnegx_plus > 1))
  expnegx = expnegx_minus*(expnegx_plus>1)+expnegx_plus*(expnegx_plus<=1)
  print(sprintf("exp(-x): %s", expnegx))
  print("tt")
  tt = -log(expnegx)
  return(tt)
}


set.seed(2023)
u_10 = runif(100)
for (alpha in 4*c(0.0625,0.125, 0.25)){
  print(sprintf("#### alpha: %s #### ", alpha))
  tt = gumbel1(u_10, alpha, 12.1)
}



N=1000
set.seed(2023)
u_N = runif(N)
y_N = rexp(N)
alpha_seq = 4*c(0.0625,0.125, 0.25)
alpha = alpha_seq[2]; print(alpha)
print(sprintf("#### alpha: %s #### ", alpha))
tt = gumbelbiexp(u_N, alpha, y_N)
Surv(tt,rep(1, 100))
s1 <- survfit(Surv(tt, rep(1,100)) ~ 1)
str(s1)
s1