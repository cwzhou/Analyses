# simulate recurrent event data with terminal event (death)

# Ghosh and Lin Section 4 (page 557): 
# “The survival time and the gap times for recurrent events were all exponential RVs. 
# To create the dependence among these variables, 
# we used Gumbel (1960)'s bivariate exponential distribution with correlation \rho. 
# To be specific, $\rho$ is the correlation between the first recurrent time and survival time 
# as well as the common correlation between any two successive gap times. 
# We set $\rho = 0,0.0625,0.125$, and $0.25$. 
# Marginally, the survival time was exponential with rate $\lambda^D = 0.25$, 
# while each gap time was exponential with rate $\lambda = 1$. 
# The censoring time was a $uniform[0,10]$ random variable, 
# which resulted in approximately two observed events per subject. 
# We set $n = 50$ and $100$. 
# For each simulation setting, $10,000$ samples were generated. 
# The summary statistics from these studies are presented in Table 1."

# https://stats.stackexchange.com/questions/318738/how-to-generate-random-samples-from-gumbel-s-bivariate-exponential-distribution

n = 50
# Simulate T marginally: T ~ exp(0.25)
T = rexp(n,0.25)
# Simulate G given this realization of T: G|T =
c1 = 
c2 = 
