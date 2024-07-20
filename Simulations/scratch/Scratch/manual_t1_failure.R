# Function definition
calculate_t1 <- function(u1, mass_p, covariate, beta_1_A) {
  print(sprintf("Zbeta1A: %s", covariate %*% beta_1_A))
  term1 <- u1 * (1-((1-mass_p)^(exp(covariate %*% beta_1_A))))
  term2 <- 1 - term1
  term3 <- log(term2) / exp(covariate %*% beta_1_A)
  t1_value <- -log(1 - mass_p^(-1) * (1 - exp(term3)))
  return(t1_value)
}


# i  = 799
u1 <- 0.77009078883566
mass_p <- 0.3
covariate <- c(1, -0.55216203,  0.04494113, -0.31196019)
action = 1
print(action)
if (action == 1){
  beta_1_A = c(0, 1 * c(3,6,5))
} else{
  beta_1_A = c(0, 1 * c(-3,5,4))
}
print(beta_1_A)

result_t <- calculate_t1(u1, mass_p, covariate, beta_1_A)
print(result_t)



# i  = 800
u1 <- 0.962314439937472
mass_p <- 0.3
covariate <- c(1, -2.328500, 1, 1)#-2.955034, -2.511810)
action = 1
print(action)
if (action == 1){
  beta_1_A = c(0, 1 * c(3,6,5))
} else{
  beta_1_A = c(0, 1 * c(-3,5,4))
}
print(beta_1_A)

result_t <- calculate_t1(u1, mass_p, covariate, beta_1_A)
print(result_t)

