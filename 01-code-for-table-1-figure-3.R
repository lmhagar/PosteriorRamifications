## code to reproduce Table 1 and Figure 3 in the main text

## load required R packages
require(ggplot2)
require(ggpubr)
require(cowplot)
require(MASS)
require(copula)
require(foreach)
require(doParallel)
require(doSNOW)

## helper function determine whether the parameter that generated the data is
## contained in the HPD set; this function also computes posterior correlation
## in terms of Kendall's tau
contains_par <- function(par, alphas, betas, prob = 0.95, seed = 1, gener = 100000, 
                         retain = 10000, rho = 0, m = 10){
  ## par: parameter value used to generate data
  ## alphas: alpha parameters for marginal beta priors for Z components
  ## betas: beta parameters for marginal beta priors for Z components
  ## prob: nominal coverage for HPD set
  ## seed: integer for reproducibility
  ## gener: number of observations to generate from proposal distribution (sampling-resampling)
  ## retain: number of observations to retain in sampling-resampling
  ## rho: Pearson's correlation for Gaussian copula in analysis prior
  ## m: sample size
  
  set.seed(seed)
  z1 <- par[1]
  z2 <- par[2]
  
  p1 <- z1
  p2 <- z2*(1-z1)
  p3 <- (1-z2)*(1-z1)
  
  ## generate multinomial data
  dat <- as.numeric(rmultinom(1, m, c(p1, p2, p3)))
  
  ## define normal copula used for the resampling weights
  norm.cop <- normalCopula(c(rho), dim = 2, dispstr = "un")
  
  ## number of initial points for sampling-resampling
  n_sim2 = gener
  
  ## using sampling-importance-resampling (Rubin, 1987) to get posterior draws at
  ## various sample sizes; the proposal distribution consists of the posterior if
  ## the marginal priors were joined with an independence copula
  v1s <- runif(n_sim2); v2s <- runif(n_sim2)
  z1s <- qbeta(v1s, dat[1] + alphas[1], dat[2] + dat[3] + betas[1])
  z2s <- qbeta(v2s, dat[2] + alphas[2], dat[3] + betas[2])
  u1s <- pbeta(z1s,alphas[1], betas[1])
  u2s <- pbeta(z2s,alphas[2], betas[2])
  
  ## convert draws from proposal distribution to the p-scale (easier for computation)
  p1s <- z1s
  p2s <- z2s*(1-z1s)
  p3s <- (1 - z2s)*(1-z1s)
  
  ## compute the numerator for the resampling weights (according to an expression that is
  ## proportional to the log-posterior) on the log-scale
  num_w <- dbeta(z1s, alphas[1], betas[1], log = TRUE) + dbeta(z2s, alphas[2], betas[2], log = TRUE) +
    dCopula(cbind(u1s, u2s), norm.cop, log = TRUE) +  
    dat[1]*log(p1s) + dat[2]*log(p2s) + dat[3]*log(p3s)
  
  ## compute the denominator of the resampling weights (according to the proposal distribution)
  ## on the log-scale
  den_w <- dbeta(z1s, dat[1] + alphas[1], dat[2] + dat[3] + betas[1], log = TRUE) + dbeta(z2s, dat[2] + alphas[2], dat[3] + betas[2], log = TRUE)
  
  ## compute resampling weights on the log-scale and exponentiate such that the weights sum to 1
  w <- num_w - den_w
  w <- w - max(w)
  w <- exp(w - max(w))
  w <- w/sum(w)
  
  ## conduct the resampling
  inds <- sample(seq(1,n_sim2,1), retain, prob = w, replace = TRUE)
  
  ## obtain the posterior draws on the Z-scale
  z1_post <- z1s[inds]
  z2_post <- z2s[inds]
  
  ## this function determines whether parameter value is in the HPD set.
  ## the HPD set is approximated using the kde2d() function in the MASS package.
  HDIbivar <- function (x, vars = 1:2, h, n = 50, lump = TRUE, prob = 0.95, 
                        xlab = NULL, ylab = NULL, lims = NULL, pt) 
  {
    parnames <- colnames(x)
    var1 <- x[, vars[1]]
    var2 <- x[, vars[2]]
    lims <- c(range(var1), range(var2))
    post1 <- kde2d(var1, var2, n = n, h = h, lims = lims)
    dx <- diff(post1$x[1:2])
    dy <- diff(post1$y[1:2])
    sz <- sort(post1$z)
    c1 <- cumsum(sz) * dx * dy
    levels <- sapply(prob, function(x) {
      approx(c1, sz, xout = 1 - x)$y
    })
    post_pt <- as.numeric(kde2d(var1, var2, n = n, h = h, lims = c(pt[1], pt[1], pt[2], pt[2]))$z)
    ifelse(post_pt >= levels, 1, 0)
  }
  
  ## this is a binary indicator determining whether the parameter value is in the HPD set
  cover <- HDIbivar(data.frame(z1 = z1_post, z2 = z2_post), n = 60, prob = prob, pt = c(z1, z2))[[1]]
  
  ## return both the binary indicator for coverage and posterior correlation
  return(c(cover, cor(z1_post, z2_post, method = "kendall")))
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

rep <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = rep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## specify parameters for marginal beta priors
alphas <- c(20,30)
betas <- c(40,30)

## specify sample sizes and correlation values
m_vals <- c(10, 100, 1000, 10000, 100000)
c_params <- round(seq(-0.95, 0.95, 0.05),2)

## generate parameter values for simulation study
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

## settings for sampling-resampling algorithm and reproducibility
gener <- 100000
retain <- 10000
dim <- 2
seed_track <- 0
for (jj in 1:length(c_params)){
  for (ii in 1:length(m_vals)){
    res_temp1 <- foreach(k=1:rep, .packages=c('copula', 'MASS'), .combine='rbind', .errorhandling = "remove",
                         .options.snow=opts) %dopar% {
                           contains_par(c(z1_sim[k], z2_sim[k]), alphas = alphas, betas = betas, seed = seed_track + k,
                                        rho = c_params[jj], m = m_vals[ii])
                         }
    seed_track <- seed_track + rep
    write.csv(res_temp1, paste0("cov_multi_", c_params[jj], "_", m_vals[ii], ".csv"), row.names = FALSE)
  }
}

## code to produce Table 1
for (ii in 1:length(m_vals)){
  tt <- read.csv(paste0("cov_multi_", rhos_keep[i], "_", m_vals[ii], ".csv"))[,2]
  round(c(min(tt), median(tt), max(tt)),4)
  
}

## put together the results into a single plot (Figure 3)
## process for m = 10
rhos_keep <- round(seq(-0.95, 0.95, 0.05),2)
samp <- 10
ctemp <- NULL
for (i in 1:length(rhos_keep)){
  ctemp[i] <- mean(read.csv(paste0("cov_multi_", rhos_keep[i], "_", samp, ".csv"))[,1])
}
c10 <- ctemp

## process for m = 100
samp <- 100
ctemp <- NULL
for (i in 1:length(rhos_keep)){
  ctemp[i] <- mean(read.csv(paste0("cov_multi_", rhos_keep[i], "_", samp, ".csv"))[,1])
}
c100 <- ctemp

## process for m = 1000
samp <- 1000
ctemp <- NULL
for (i in 1:length(rhos_keep)){
  ctemp[i] <- mean(read.csv(paste0("cov_multi_", rhos_keep[i], "_", samp, ".csv"))[,1])
}
c1000 <- ctemp

## process for m = 10000
samp <- 10000
ctemp <- NULL
for (i in 1:length(rhos_keep)){
  ctemp[i] <- mean(read.csv(paste0("cov_multi_", rhos_keep[i], "_", samp, ".csv"))[,1])
}
c10000 <- ctemp

## process for m = 100000
samp <- 100000
ctemp <- NULL
for (i in 1:length(rhos_keep)){
  ctemp[i] <- mean(read.csv(paste0("cov_multi_", rhos_keep[i], "_", samp, ".csv"))[,1])
}
c100000 <- ctemp

## create data frame with all sample sizes
df_cov <- data.frame(rho = rep(rhos_keep, 5), y = c(c10, c100, c1000, c10000, c100000),
                     n = rep(c(10, 100, 1000, 10000, 100000), each = length(rhos_keep)))

## use colour-blind-friendly palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## code to produce Figure 3
fig_cov <- ggplot(df_cov, aes(x=rho)) + theme_bw() +
  geom_line(aes(y = y, color=as.factor(n), linetype = as.factor(n)), size = 1.05) +
  labs(x= bquote('Analysis Prior '*rho), y= bquote('Coverage')) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  theme(legend.position="bottom") +
  scale_color_manual(name = bquote(italic(n)),
                     labels = c(bquote(10^1), bquote(10^2), bquote(10^3), 
                                bquote(10^4), bquote(10^5)),
                     values = c(cbb[2], cbb[3], cbb[4], cbb[7], cbb[8])) +
  scale_linetype_manual(name = bquote(italic(n)),
                        labels = c(bquote(10^1), bquote(10^2), bquote(10^3), 
                                   bquote(10^4), bquote(10^5)),
                        values = c(1,2,4,5,6)) +
  theme(legend.text=element_text(size=16), legend.title=element_text(size=16)) +
  theme(legend.key.width = unit(0.9, 'cm')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.95, lty = 3) + geom_vline(xintercept = -0.9, lty = 3)

## output as .pdf file for the article
pdf(file = "FigureCov.pdf",   # The directory you want to save the file in
    width = 5.5, # The width of the plot in inches (12.41)
    height = 4.5) # The height of the plot in inches (10.7)

fig_cov

dev.off()