library(rjags)
library(coda)
library(tidyverse)
library(gridExtra)
library(latex2exp)

flu <- read.table("flu.txt", header = TRUE)

ncountry <- length(unique(flu$Country))
ntest <- nrow(flu)/ncountry

ezk <- matrix(flu$EZK, nrow = ntest)
infected <- matrix(flu$Infected, nrow = ntest)

mydata <- list(
  ncountry = ncountry,
  ntest = ntest, 
  ezk = ezk, 
  inf = infected
)

# initial values
myinit <- list(
  alpha = rep(0, ncountry), 
  beta = rep(0, ncountry), 
  mu_a = 0, 
  mu_b = 0, 
  pres_sigma2_a = 1, 
  pres_sigma2_b = 1
)

# iterations
outiters <- 10000
nadapt <- 5000
nchains <- 1
nburn <- 5000
niters <- outiters + nburn

# Specify JAGS model:
mod = "model {

  for (j in 1:ncountry) {
    for (i in 1:ntest) {
      inf[i,j] ~ dbern(theta[i,j])
      
      logit(theta[i,j]) = alpha[j] + beta[j] * ezk[i,j]
    }
    
   alpha[j] ~ dnorm(mu_a, pres_sigma2_a)
   beta[j] ~ dnorm(mu_b, pres_sigma2_b)
    
  }
  
  mu_a ~ dnorm(0, 1/9)
  mu_b ~ dnorm(0, 1/9)
  pres_sigma2_a ~ dgamma(1, 10)
  pres_sigma2_b ~ dgamma(1, 10)
  
  
  # for tracing
  sigma2_a = 1/pres_sigma2_a
  sigma2_b = 1/pres_sigma2_b
  
}"

fit <- jags.model(
  textConnection(mod), 
  data=mydata, 
  inits=myinit, 
  n.chains=nchains, 
  n.adapt=nadapt, 
  quiet = TRUE
  )

fit.samples <- coda.samples(
  fit,
  c("mu_a", "mu_b", "sigma2_a", "sigma2_b", "alpha", "beta"), 
  n.iter=niters
  )

samples <- data.frame(fit.samples[[1]][-(1:nburn),])
# for (i in 1:nchains) {
#   samples <- rbind(
#     samples,
#     data.frame(fit.samples[[i]][-(1:nburn),], chain = i)
#   )
# }

############### figures #############

table_1 <- matrix(
  c(
    0, 0, 1, 1, 0, 1, 0, 1,
    flu %>% 
      filter(Infected == 0, EZK == 0) %>% 
      nrow(),
    flu %>% 
      filter(Infected == 0, EZK == 1) %>% 
      nrow(),
    flu %>% 
      filter(Infected == 1, EZK == 0) %>% 
      nrow(),
    flu %>% 
      filter(Infected == 1, EZK == 1) %>% 
      nrow()
  ),
  nrow = 4,
  dimnames = list(NULL , c("Infected","EZK", "N"))
)


fig_1.1 <- samples %>% 
  pivot_longer(cols = colnames(samples)[1:10], names_to = "ja", values_to = "alpha") %>% 
  pivot_longer(cols = colnames(samples)[11:20], names_to = "jb", values_to = "beta") %>% 
  ggplot() +
    geom_density(aes(x = alpha, color = ja), show.legend = F) +
    geom_density(aes(x = beta, color = jb), show.legend = F) +
    #geom_density(aes(x = mu_a), linetype = "dashed") +
    #geom_density(aes(x = mu_b), linetype = "dashed") +
    xlab(TeX("$\\alpha_j$ and $\\beta_j$ for $j=1,...,10$")) +
    theme_minimal()

fig_1.3 <- samples %>% 
  ggplot() +
    geom_density(aes(x = sigma2_a), color = "red") +
    geom_density(aes(x = sigma2_b), color = "blue") +
    xlab(TeX("$\\sigma^2_\\alpha$ and $\\sigma^2_\\beta$")) +
    theme_minimal() +
    labs(caption = "Figure 1")

fig_1.2 <- samples %>% 
  pivot_longer(cols = colnames(samples)[1:10], names_to = "ja", values_to = "alpha") %>% 
  pivot_longer(cols = colnames(samples)[11:20], names_to = "jb", values_to = "beta") %>% 
  select(alpha, beta, ja, jb) %>% 
  mutate(iteration = row_number()/100) %>% 
  ggplot() +
    geom_line(aes(x = iteration, y = alpha, color = ja), show.legend = F) +
    geom_line(aes(x = iteration, y = beta, color = jb), show.legend = F) +
    ylab(TeX("$\\alpha_j$ and $\\beta_j$")) +
    theme_minimal()

means <- colMeans(samples)
table_2 <- tibble(
  country = LETTERS[1:10],
  alpha.mean = means[1:10],
  beta.mean = means[11:20],
  logit.theta.EZKpos = alpha.mean + beta.mean,
  logit.theta.EZKneg = alpha.mean,
  infect.odds.EZKpos = exp(logit.theta.EZKpos),
  infect.odds.EZKneg = exp(logit.theta.EZKneg),
  infect.prob.EZKpos = infect.odds.EZKpos / (1 + infect.odds.EZKpos),
  infect.prob.EZKneg = infect.odds.EZKneg / (1 + infect.odds.EZKneg)
) %>% 
  select(country, alpha.mean, beta.mean, infect.prob.EZKpos, infect.prob.EZKneg)

fig_2 <- table_2 %>% 
  mutate(`false negative` = infect.prob.EZKneg, `false positive` = 1 - infect.prob.EZKpos) %>% 
  select(`false negative`, `false positive`, country) %>% 
  pivot_longer(cols = c("false negative", "false positive"), names_to = "event", values_to = "probability") %>% 
  ggplot(aes(x = country, y = probability, fill = event)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    labs(caption = "Figure 2", title = "Probability of False Results")

















