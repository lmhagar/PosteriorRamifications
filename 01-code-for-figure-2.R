## code to reproduce Figure 2 in the main text

## load required R packages
require(ggplot2)
require(ggpubr)
require(cowplot)
require(MASS)
require(copula)

## give the conditional median probabilities on Z-scale (theta_0 values)
medians <- c(1/3, 1/2)

## define parameters for marginal beta priors
alphas <- c(20,40)
betas <- c(30,30)

## define correlation matrix for Gaussian copula with rho = -0.9
c_param <- -0.9
cop_mat <- matrix(c_param, nrow = length(medians), ncol = length(medians)) + (1 - c_param)*diag(length(medians))

## simulate 10000 observations from Gaussian copula
n_sim = 10000
set.seed(2025)
obs <- pnorm(mvrnorm(n = n_sim, mu = rep(0, length(medians)), Sigma = cop_mat))

## convert to 10000 observations from the joint prior on Z1 and Z2
z1_sim <- qbeta(obs[,1], alphas[1], betas[1])
z2_sim <- qbeta(obs[,2], alphas[2], betas[2])

## compute the induced prior draws on p1, p2, and p3
p1_sim <- z1_sim
p2_sim <- z2_sim*(1-z1_sim)
p3_sim <- (1-z2_sim)*(1-z1_sim)

## simulate 10000 observations according to these prior draws for the multinomial probabilities
set.seed(2)
dat <- NULL
for (i in 1:length(p1_sim)){
  dat <- rbind(dat, as.numeric(rmultinom(1, 1, cbind(p1_sim[i], p2_sim[i], p3_sim[i]))))
}

## store the data to make the figures at different sample sizes
dat_store <- dat

## using sampling-importance-resampling (Rubin, 1987) to get posterior draws at
## various sample sizes; the proposal distribution consists of the two marginal priors
## joined by the independence copula
set.seed(3)
n_sim2 = 1000000
u1s <- runif(n_sim2); u2s <- runif(n_sim2)
z1s <- qbeta(u1s,alphas[1], betas[1])
z2s <- qbeta(u2s,alphas[2], betas[2])

## convert draws from proposal distribution to the p-scale (easier for computation)
p1s <- z1s
p2s <- z2s*(1-z1s)
p3s <- (1 - z2s)*(1-z1s)

## define normal copula used for the resampling weights
norm.cop <- normalCopula(c(-0.9), dim = 2, dispstr = "un")

samps <- c(10, 100, 1000, 10000)
for (k in 1:length(samps)){
  ## take the appropriate number of observations
  dat <- colSums(dat_store[1:samps[k],])
  
  ## compute the numerator for the resampling weights (according to an expression that is
  ## proportional to the log-posterior) on the log-scale
  num_w <- dbeta(z1s, alphas[1], betas[1], log = TRUE) + dbeta(z2s, alphas[2], betas[2], log = TRUE) +
    dCopula(cbind(u1s, u2s), norm.cop, log = TRUE) +  
    dat[1]*log(p1s) + dat[2]*log(p2s) + dat[3]*log(p3s)
  
  ## compute the denominator of the resampling weights (according to the proposal distribution)
  ## on the log-scale
  den_w <- dbeta(z1s, alphas[1], betas[1], log = TRUE) + dbeta(z2s, alphas[2], betas[2], log = TRUE)
  
  ## compute resampling weights on the log-scale and exponentiate such that the weights sum to 1
  w <- num_w - den_w
  w <- w - max(w)
  w <- exp(w - max(w))
  w <- w/sum(w)
  
  ## conduct the resampling
  set.seed(k + 3)
  inds <- sample(seq(1,1000000,1), 10000, prob = w, replace = TRUE)
  
  ## obtain the posterior draws on the Z-scale
  z1_post <- z1s[inds]
  z2_post <- z2s[inds]
  ## output the posterior correlation in terms of Kendall's tau
  print(round(cor(z1_post, z2_post, method = "kendall"),3))
  
  ## save posterior draws for each sample size
  assign(paste0("z1_post", samps[k]), z1_post)
  assign(paste0("z2_post", samps[k]), z2_post)
}

## estimate copula for prior using ranks
dat_prior <- data.frame(u1 = (rank(z1_sim)-0.5)/10000, u2 = (rank(z2_sim)-0.5)/10000)
plot1 <- ggplot(dat_prior, aes(x=u1, y=u2)) + theme_bw() +
  geom_point(size=2, alpha = 0.12) +
  labs(title = bquote("Prior:"~tau~"="~-0.713)) +
  labs(x=bquote('\n'~"F"[1]*"("*italic(Z)[1]*")"), y=bquote("F"[2]*"("*italic(Z)[2]*")"~'\n')) +
  theme(plot.title = element_text(hjust = 0.5,size=16,margin=unit(c(0,0,5,0), "mm"))) +
  theme(plot.subtitle = element_text(hjust = 0.5,size=16)) +
  theme(axis.text=element_text(size=14),
        axis.title.x=element_text(size=16, margin=unit(c(3,0,0,0), "mm")),
        axis.title.y=element_text(size = 16, margin=unit(c(0,3,0,0), "mm")))

## estimate copula for posterior using ranks; repeat for n = 100, 1000, and 10000
dat_10 <- data.frame(u1 = (rank(z1_post10)-0.5)/10000, u2 = (rank(z2_post10)-0.5)/10000)
plot2 <- ggplot(dat_10, aes(x=u1, y=u2)) + theme_bw() +
  geom_point(size=2, alpha = 0.12) +
  labs(title = bquote("Posterior:"~italic(n)~"="~10^1*","~tau~"="~-0.683)) +
  labs(x=bquote('\n'~"F"[1]*"("*italic(Z)[1]*")"), y=bquote("F"[2]*"("*italic(Z)[2]*")"~'\n')) +
  theme(plot.title = element_text(hjust = 0.5,size=16,margin=unit(c(0,0,5,0), "mm"))) +
  theme(plot.subtitle = element_text(hjust = 0.5,size=16)) +
  theme(axis.text=element_text(size=14),
        axis.title.x=element_text(size=16, margin=unit(c(3,0,0,0), "mm")),
        axis.title.y=element_text(size = 16, margin=unit(c(0,3,0,0), "mm")))

dat_100 <- data.frame(u1 = (rank(z1_post100)-0.5)/10000, u2 = (rank(z2_post100)-0.5)/10000)
plot3 <- ggplot(dat_100, aes(x=u1, y=u2)) + theme_bw() +
  geom_point(size=2, alpha = 0.12) +
  labs(title = bquote("Posterior:"~italic(n)~"="~10^2*","~tau~"="~-0.501)) +
  labs(x=bquote('\n'~"F"[1]*"("*italic(Z)[1]*")"), y=bquote("F"[2]*"("*italic(Z)[2]*")"~'\n')) +
  theme(plot.title = element_text(hjust = 0.5,size=16,margin=unit(c(0,0,5,0), "mm"))) +
  theme(plot.subtitle = element_text(hjust = 0.5,size=16)) +
  theme(axis.text=element_text(size=14),
        axis.title.x=element_text(size=16, margin=unit(c(3,0,0,0), "mm")),
        axis.title.y=element_text(size = 16, margin=unit(c(0,3,0,0), "mm")))

dat_1000 <- data.frame(u1 = (rank(z1_post1000)-0.5)/10000, u2 = (rank(z2_post1000)-0.5)/10000)
plot4 <- ggplot(dat_1000, aes(x=u1, y=u2)) + theme_bw() +
  geom_point(size=2, alpha = 0.12) +
  labs(title = bquote("Posterior:"~italic(n)~"="~10^3*","~tau~"="~-0.168)) +
  labs(x=bquote('\n'~"F"[1]*"("*italic(Z)[1]*")"), y=bquote("F"[2]*"("*italic(Z)[2]*")"~'\n')) +
  theme(plot.title = element_text(hjust = 0.5,size=16,margin=unit(c(0,0,5,0), "mm"))) +
  theme(plot.subtitle = element_text(hjust = 0.5,size=16)) +
  theme(axis.text=element_text(size=14),
        axis.title.x=element_text(size=16, margin=unit(c(3,0,0,0), "mm")),
        axis.title.y=element_text(size = 16, margin=unit(c(0,3,0,0), "mm")))

dat_10000 <- data.frame(u1 = (rank(z1_post10000)-0.5)/10000, u2 = (rank(z2_post10000)-0.5)/10000)
plot5 <- ggplot(dat_10000, aes(x=u1, y=u2)) + theme_bw() +
  geom_point(size=2, alpha = 0.12) +
  labs(title = bquote("Posterior:"~italic(n)~"="~10^4*","~tau~"="~-0.018)) +
  labs(x=bquote('\n'~"F"[1]*"("*italic(Z)[1]*")"), y=bquote("F"[2]*"("*italic(Z)[2]*")"~'\n')) +
  theme(plot.title = element_text(hjust = 0.5,size=16,margin=unit(c(0,0,5,0), "mm"))) +
  theme(plot.subtitle = element_text(hjust = 0.5,size=16)) +
  theme(axis.text=element_text(size=14),
        axis.title.x=element_text(size=16, margin=unit(c(3,0,0,0), "mm")),
        axis.title.y=element_text(size = 16, margin=unit(c(0,3,0,0), "mm")))

## combine subplots into grid of two rows
fig.row1 <- plot_grid(plot1 + theme(plot.margin=unit(c(0.45,0.45,0.45,0.45),"cm")), 
                      plot2 + theme(plot.margin=unit(c(0.45,0.45,0.45,0.45),"cm")),
                      plot3 + theme(plot.margin=unit(c(0.45,0.45,0.45,0.45),"cm")),
                      rel_widths = c(1, 1,1), ncol = 3)
fig.row2 <- plot_grid(plot4 + theme(plot.margin=unit(c(0.45,0.45,0.45,0.45),"cm")), 
                      plot5 + theme(plot.margin=unit(c(0.45,0.45,0.45,0.45),"cm")),
                      NULL,
                      rel_widths = c(1, 1,1), ncol = 3)
fig <- plot_grid(fig.row1, fig.row2, nrow = 2)

pdf(file = "Fig2_PR.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches (12.41)
    height = 8) # The height of the plot in inches (10.7)

fig
dev.off()