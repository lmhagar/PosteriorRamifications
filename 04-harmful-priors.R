## code to reproduce analysis with harmful priors for reviewer response

## load required R packages
require(ggplot2)
require(ggpubr)
require(cowplot)
require(MASS)
require(copula)
require(foreach)
require(doParallel)
require(doSNOW)

## define parameters for marginal beta priors
alphas <- c(20,30)
betas <- c(40,30)

## get parameter values for simulation study
dim <- 2
## define correlation matrix for Gaussian copula with rho = -0.9
c_param <- -0.9
cop_mat <- matrix(c_param, nrow = dim, ncol = dim) + (1 - c_param)*diag(dim)

## simulate 10000 observations from Gaussian copula
n_sim = 10000
set.seed(2025)
obs <- pnorm(mvrnorm(n = n_sim, mu = rep(0, dim), Sigma = cop_mat))

## convert to 10000 observations from the joint prior on Z1 and Z2
z1_sim <- qbeta(obs[,1], alphas[1], betas[1])
z2_sim <- qbeta(obs[,2], alphas[2], betas[2])

## now get the prior parameters to put into the for loop
rhos <- seq(-0.95, 0.95, 0.05)
Sigma_vec <- NULL
mu_vec <- NULL
Rs <- NULL
for (j in 1:length(rhos)){
  alphas <- c(20,40)
  betas <- c(30,30)
  
  dim <- 2
  ## define correlation matrix for Gaussian copula with rho = -0.9
  c_param <- rhos[j]
  cop_mat <- matrix(c_param, nrow = dim, ncol = dim) + (1 - c_param)*diag(dim)
  
  ## simulate 10000 observations from Gaussian copula
  n_sim = 200000
  set.seed(2025)
  obs <- pnorm(mvrnorm(n = n_sim, mu = rep(0, dim), Sigma = cop_mat))
  
  ## convert to 10000 observations from the joint prior on Z1 and Z2
  z1_sim <- qbeta(obs[,1], alphas[1], betas[1])
  z2_sim <- qbeta(obs[,2], alphas[2], betas[2])
  
  ## approximate priors with normal distributions using simulation
  mu_pi <- c(mean(z1_sim), mean(z2_sim))
  Sigma_pi <- c(var(z1_sim), cov(z1_sim, z2_sim), cov(z1_sim, z2_sim), var(z2_sim))
  
  ## mu_vec contains mu_pi for each prior setting from equation 11 in Reimherr et al.
  mu_vec <- rbind(mu_vec, mu_pi)
  
  ## Sigma_vec contains sigma_pi for each prior setting from equation 11 in Reimherr et al.
  Sigma_vec <- rbind(Sigma_vec, Sigma_pi)
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

rep <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = rep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## repeat simulation settings for Figure 3
alphas <- c(20,30)
betas <- c(40,30)
m_vals <-  c(10, 100, 1000, 10000, 100000)
c_param <- -0.9
gener <- 100000
retain <- 10000
dim <- 2

## iterate over each sample size to create the supplementary plot
for (ii in 1:length(m_vals)){
  res_temp1 <- foreach(k=1:rep, .packages=c('copula', 'MASS'), .combine='rbind',
                       .options.snow=opts) %dopar% {
                         
                         ## generate a sample
                         set.seed(rep + 10000*(ii - 1))
                         z1 <- z1_sim[k]
                         z2 <- z2_sim[k]
                         m <- m_vals[ii]
                         
                         p1 <- z1
                         p2 <- z2*(1-z1)
                         p3 <- (1-z2)*(1-z1)
                         
                         dat <- as.numeric(rmultinom(1, m_val, c(p1, p2, p3)))
                         
                         ## join the marginal priors with an independence copula for the baseline prior
                         c_prior <- 0
                         ## define normal copula used for the resampling weights
                         norm.cop <- normalCopula(c(c_prior), dim = 2, dispstr = "un")
                         
                         ## number of points to sample (we do not need sampling-resampling when using 
                         ## the independence copula)
                         n_sim2 <- retain
                         
                         ## use uninformative GAMMA(2,1) priors as the baseline prior from
                         ## equation 11 in Reimherr et al.
                         alphas <- c(1,1)
                         betas <- c(2,1)
                         
                         ## using sampling-importance-resampling (Rubin, 1987) to get posterior draws at
                         ## various sample sizes; the proposal distribution consists of the posterior if
                         ## the marginal priors were joined with an independence copula
                         v1s <- runif(n_sim2); v2s <- runif(n_sim2)
                         z1s <- qbeta(v1s, dat[1] + alphas[1], dat[2] + dat[3] + betas[1])
                         z2s <- qbeta(v2s, dat[2] + alphas[2], dat[3] + betas[2]
                         
                         ## do resampling is required
                         z1_post <- z1s
                         z2_post <- z2s
                         
                         mu_hat <- c(mean(z1_post), mean(z2_post))
                         
                         Sigma_hat <- rbind(c(var(z1_post), cov(z1_post, z2_post)),
                                            c(cov(z1_post, z2_post), var(z2_post)))
                         
                         ## compute the relative information/improvement factor (R(n)) for 
                         ## all prior dependence structures using the same baseline posterior sample
                         rhos <- seq(-0.95, 0.95, 0.05)
                         Rs <- NULL
                         for (jj in 1:length(rhos)){
                           
                           ## extract the proper mu_pi and Sigma_pi values for each setting
                           mu_pi <- c(mu_vec[jj,1], mu_vec[jj,2])
                           
                           Sigma_pi <- rbind(c(Sigma_vec[jj,1], Sigma_vec[jj,2]),
                                             c(Sigma_vec[jj,3], Sigma_vec[jj,4]))
                           
                           ## compute r_hat from Reimherr et al. (d = 2)
                           r_hat <- 0.5*sum(diag((Sigma_hat%*%solve(Sigma_pi))))
                           
                           ## compute Delta_hat from Reimherr et al.
                           Delta_hat <- as.numeric(t((mu_hat - mu_pi))%*%solve(Sigma_pi)%*%(mu_hat - mu_pi))
                           
                           ## compute R_hat from Reimherr et al.
                           R_hat <- r_hat*(1 - (2/(Delta_hat - 1) + r_hat/(1 + r_hat))^(-1))
                           Rs[jj] <- R_hat
                         }
                         Rs
                       }
  
  write.csv(res_temp1, paste0("harmful_", m_vals[ii], ".csv"), row.names = FALSE)
}

## create the composite figure
# output as .pdf file for the article
pdf(file = "FigureHarmful.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches (12.41)
    height = 8) # The height of the plot in inches (10.7)

par(mfrow = c(3, 2), mar = c(4, 4.5, 2, 0.5))

h10 = read.csv("harmful_10.csv")
plot(rhos, rep(0, length(rhos)), col = "white", type = "l",
     ylab = "Relative Improvement",
     xlab = bquote("Analysis Prior "*rho*" ("*italic(n)*" = "*10^1*")"), ylim = c(-2, 12))
for (i in 1:1000){
  lines(rhos, h10[i,], col = adjustcolor("grey", 0.1))
}
lines(rhos, apply(h10, 2, median), type = "l", col = "firebrick")

h100 = read.csv("harmful_100.csv")
plot(rhos, rep(0, length(rhos)), col = "white", type = "l",
     ylab = "Relative Improvement",
     xlab = bquote("Analysis Prior "*rho*" ("*italic(n)*" = "*10^2*")"),
     ylim = c(-2, 12))
for (i in 1:1000){
  lines(rhos, h100[i,], col = adjustcolor("grey", 0.1))
}
lines(rhos, apply(h100, 2, median), type = "l", col = "firebrick")

h1000 = read.csv("harmful_1000.csv")
plot(rhos, rep(0, length(rhos)), col = "white", type = "l",
     ylab = "Relative Improvement",
     xlab = bquote("Analysis Prior "*rho*" ("*italic(n)*" = "*10^3*")"), ylim = c(-1, 2))
for (i in 1:1000){
  lines(rhos, h1000[i,], col = adjustcolor("grey", 0.1))
}
lines(rhos, apply(h1000, 2, median), type = "l", col = "firebrick")

h10000 = read.csv("harmful_10000.csv")
plot(rhos, rep(0, length(rhos)), col = "white", type = "l",
     ylab = "Relative Improvement",
     xlab = bquote("Analysis Prior "*rho*" ("*italic(n)*" = "*10^4*")"), ylim = c(-1, 0.5))
for (i in 1:1000){
  lines(rhos, h10000[i,], col = adjustcolor("grey", 0.1))
}
lines(rhos, apply(h10000, 2, median), type = "l", col = "firebrick")

h100000 = read.csv("harmful_100000.csv")
plot(rhos, rep(0, length(rhos)), col = "white", type = "l",
     ylab = "Relative Improvement",
     xlab = bquote("Analysis Prior "*rho*" ("*italic(n)*" = "*10^5*")"), ylim = c(-0.75, 0.1))
for (i in 1:1000){
  lines(rhos, h100000[i,], col = adjustcolor("grey", 0.1))
}
lines(rhos, apply(h100000, 2, median), type = "l", col = "firebrick")

par(mfrow = c(1, 1))

dev.off()