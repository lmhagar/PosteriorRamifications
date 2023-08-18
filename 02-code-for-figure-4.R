## code to reproduce Figure 4 in the main text

## load required R packages
require(foreach)
require(doParallel)
require(doSNOW)
require(nleqslv)
require(MASS)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(scales)

## set up parallelization for m = 10000 repetitions per
## sample size and scenario combination
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## helper function to find the modes for the regression model
## x is the theta = (beta_1, beta_2) combination
## cop = TRUE for t-copula and FALSE for independence copula
## dat is a data frame with columns for y (response) and covariates x1 and x2
fn <- function(x, cop = TRUE, dat) {
  
  yy <- dat$y
  xx1 <- dat$x1
  xx2 <- dat$x2
  
  ## t-copula has 4 degrees of freedom
  nu <- 4
  zz <- x
  
  ## prior variance for the regression coefficients is 1
  ## variance for the error term of the response is 5
  varb <- 1
  vary <- 5
  
  ## convert from the theta scale to the (0,1)-scale
  ## if u in (0,1) > 0.5, we take 1-u. This is more stable because
  ## 1 - u may round to 1. This is relevant for case 6. This is also 
  ## appropriate because the t-copula with diagonal correlation matrix
  ## and the independence copula are symmetric
  uu <- pnorm(zz[1], 0, sqrt(varb))
  if (uu > 0.5){
    uu <- pnorm(zz[1], 0, sqrt(varb), lower.tail = FALSE)
  }
  vv <- pnorm(zz[2], 0, sqrt(varb))
  if (vv > 0.5){
    vv <- pnorm(zz[2], 0, sqrt(varb), lower.tail = FALSE)
  }
  
  if (cop == TRUE){
    ## first term is the derivative with respect to theta1 of the logarithm of the marginal prior
    dlogpdz11 <- -1*zz[1]/varb 
    
    ## second two terms comprise the derivative with respect to theta1 of the logarithm of the t-copula density
    ## this is 0 for the independence copula
    dlogpdz12 <- -0.5*(nu + 2)*(1 + qt(uu, nu)^2/nu + qt(vv, nu)^2/nu)^(-1)*((2/nu)*qt(uu, nu))*(1/dt(qt(uu, nu), nu))*dnorm(zz[1], 0, sqrt(varb))
    dlogpdz13 <- 0.5*(nu + 1)*(1 + qt(uu, nu)^2/nu)^(-1)*((2/nu)*qt(uu, nu))*(1/dt(qt(uu, nu), nu))*dnorm(zz[1], 0, sqrt(varb))
    
    ## combine the derivative of the additive functions on the log scale
    dlogpdz1 <- dlogpdz11 + dlogpdz12 + dlogpdz13
    
    ## we need to instead multiply the last two components by negative 1 if we took 1-u
    if (pnorm(zz[1], 0, sqrt(varb)) > 0.5){
      dlogpdz1 <- dlogpdz11 - dlogpdz12 - dlogpdz13
    }
  }
  else{
    ## do not take derivatives of last two components for the independence copula
    dlogpdz1 <- -1*zz[1]/varb
  }
  
  ## compute derivatives of log-likelihood with respect to theta1
  dlogLdz1 <- (sum(yy*xx1) - zz[2]*sum(xx1*xx2) - zz[1]*sum(xx1^2))/vary
  
  ## combine derivatives of log-likelihood with logarithm of prior density
  dz1 <- dlogpdz1 + dlogLdz1
  
  ## repeat this process for theta2
  if (cop == TRUE){
    dlogpdz21 <- -1*zz[2]/varb 
    dlogpdz22 <- -0.5*(nu + 2)*(1 + qt(uu, nu)^2/nu + qt(vv, nu)^2/nu)^(-1)*((2/nu)*qt(vv, nu))*(1/dt(qt(vv, nu), nu))*dnorm(zz[2], 0, sqrt(varb))
    dlogpdz23 <- 0.5*(nu + 1)*(1 + qt(vv, nu)^2/nu)^(-1)*((2/nu)*qt(vv, nu))*(1/dt(qt(vv, nu), nu))*dnorm(zz[2], 0, sqrt(varb))
    
    dlogpdz2 <- dlogpdz21 + dlogpdz22 + dlogpdz23
    
    if (pnorm(zz[2], 0, sqrt(varb)) > 0.5){
      dlogpdz2 <- dlogpdz21 - dlogpdz22 - dlogpdz23
    }
    
  }
  else{
    dlogpdz2 <- -1*zz[2]/varb
  }
  dlogLdz2 <- (sum(yy*xx2) - zz[1]*sum(xx1*xx2) - zz[2]*sum(xx2^2))/vary
  
  dz2 <- dlogpdz2 + dlogLdz2
  
  ## return vector of partial derivatives of logarithm of posterior with respect to theta1 and theta2
  ## these should be 0 at the posterior mode
  return(c(dz1, dz2))
}

## set (u1, u2) values for first five settings
u1 <- c(0.5, 0.495, pt(1,4), 0.5, 0.85)
u2 <- c(0.5, 0.495, pt(1,4), 0.99, 0.9)

## convert to coefficient (beta) scale
q1 <- qnorm(u1, 0, 1)
q2 <- qnorm(u2, 0, 1)

## add coefficients for final setting (case 6)
q1 <- c(q1, -5)
q2 <- c(q2, 8)

## define vector of sample sizes for simulation
samps <- c(5,10,15,20,25, 50, 75, 100, 125, 150, 175, 200, seq(250, 950,50), 
           seq(1000,10000,100), seq(10000,25000,500), seq(30000,100000,5000))

## repeat for each of the six scenarios
for (i in 1:length(q1)){
  res <- NULL
  for (j in 1:length(samps)){
    ## repeat process 10000 times in parallel for each sample size and scenario combination
    set.seed(10000*j, kind = "L'Ecuyer-CMRG")
    res_temp <- foreach(k=1:m, .packages=c('nleqslv', 'MASS'), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          ## simulate covariates x1 and x2
                          temp = mvrnorm(n=samps[j], mu = c(0,0), Sigma = rbind(c(1, 0), c(0,1)))
                          df = data.frame(x1 = temp[,1], x2 = temp[,2])
                          
                          ## simulate response y and create data frame
                          y = q1[i]*df$x1 + q2[i]*df$x2 + rnorm(samps[j],0, sqrt(5))
                          dat_temp = cbind(df, y)
                          
                          ## compute the posterior mode with each copula for the same sample
                          modet <- nleqslv(c(q1[i], q2[i]), fn, cop = TRUE, dat = dat_temp)$x
                          modeind <- nleqslv(c(q1[i], q2[i]), fn, cop = FALSE, dat = dat_temp)$x
                          
                          ## compute the squared Euclidean distance (L2-norm) between each mode and the
                          ## true coefficients
                          dist <- (modet[1] - q1[i])^2 + (modet[2] - q2[i])^2
                          disind <- (modeind[1] - q1[i])^2 + (modeind[2] - q2[i])^2
                          
                          ## return modes and distances from true coefficients
                          c(modet, dist, modeind, disind)
                        }
    ## for each sample size, estimate the probability of the t-copula's posterior mode
    ## being closer to theta_0 than the mode for the independence copula. Also estimate the
    ## expected absolute difference between the L2-norms (D2 - D1) in the text
    res <- rbind(res, c(samps[j], mean(ifelse(res_temp[,3] < res_temp[,6],1,0)), 
                        mean(abs(sqrt(res_temp[,3]) - sqrt(res_temp[,6])))))
    
    ## output the results to a .csv file for each scenario
    write.csv(res, paste0("reg_modes_", i, ".csv"), row.names = FALSE)
  }
}

## plots for Figure 4

## start with scenario 1 and define scale for right vertical axis
df1 <- read.csv("reg_modes_1.csv")
scale = max(df1$V3)

## create plot; ensure the horizontal axis is on the logarithmic scale
plot1 <- ggplot(df1, aes(x = V1, y = V2)) + theme_bw() +
  geom_line(aes(color = "col1", linetype = "col1"), lwd = 1.25) +
  geom_line(aes(y = V3/scale, color = "col2", linetype = "col2"), lwd = 1.25) +
  labs(title = bquote("Case 1: ("*italic(u)["1,0"]*","*italic(u)["2,0"]*")"~"= ("*0.5*","*0.5*")")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name=bquote('|'*italic(D)[2] - italic(D)[1]*'|'))) +
  labs(x=bquote('\n'~italic(n)), 
       y=bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')')) +
  scale_color_manual(values = c("firebrick4", "steelblue4"), name="",
                     labels = c("Probability", "Mean Difference")) +
  scale_linetype_manual(values = c(1,2), name="",
                        labels = c("Probability", "Mean Difference")) +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(size = 17, margin=margin(10,0,0,0))) +
  theme(axis.title.y = element_text(size = 17, margin=margin(0,10,0,0), colour = "firebrick4")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15, colour = "firebrick4")) +
  theme(plot.title=element_text(size=18)) +
  theme(plot.margin=unit(c(5.5,18.5,5.5,5.5), "pt")) +
  theme(axis.title.y.right = element_text(vjust=3, colour = "steelblue4")) +
  theme(axis.text.y.right = element_text(colour = "steelblue4"))

## repeat process for case 2
df2 <- read.csv("reg_modes_2.csv")
scale2 = max(df2$V3)

plot2 <- ggplot(df2, aes(x = V1, y = V2)) + theme_bw() +
  geom_line(aes(color = "col1", linetype = "col1"), lwd = 1.25) +
  geom_line(aes(y = V3/scale2, color = "col2", linetype = "col2"), lwd = 1.25) +
  labs(title = bquote("Case 2: ("*italic(u)["1,0"]*","*italic(u)["2,0"]*")"~"= ("*0.495*","*0.495*")")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(sec.axis = sec_axis(~.*scale2, name=bquote('|'*italic(D)[2] - italic(D)[1]*'|'))) +
  labs(x=bquote('\n'~italic(n)), 
       y=bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')')) +
  scale_color_manual(values = c("firebrick4", "steelblue4"), name="",
                     labels = c("Probability", "Mean Difference")) +
  scale_linetype_manual(values = c(1,2), name="",
                        labels = c("Probability", "Mean Difference")) +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(size = 17, margin=margin(10,0,0,0))) +
  theme(axis.title.y = element_text(size = 17, margin=margin(0,10,0,0), colour = "firebrick4")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15, colour = "firebrick4")) +
  theme(plot.title=element_text(size=18)) +
  theme(plot.margin=unit(c(5.5,18.5,5.5,5.5), "pt")) +
  theme(axis.title.y.right = element_text(vjust=3, colour = "steelblue4")) +
  theme(axis.text.y.right = element_text(colour = "steelblue4"))

## repeat process for case 3
df3 <- read.csv("reg_modes_3.csv")
scale3 = max(df3$V3)

plot3 <- ggplot(df3, aes(x = V1, y = V2)) + theme_bw() +
  geom_line(aes(color = "col1", linetype = "col1"), lwd = 1.25) +
  geom_line(aes(y = V3/scale3, color = "col2", linetype = "col2"), lwd = 1.25) +
  labs(title = bquote("Case 3: ("*italic(u)["1,0"]*","*italic(u)["2,0"]*")"~"= ("*0.813*","*0.813*")")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(sec.axis = sec_axis(~.*scale3, name=bquote('|'*italic(D)[2] - italic(D)[1]*'|'))) +
  labs(x=bquote('\n'~italic(n)), 
       y=bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')')) +
  scale_color_manual(values = c("firebrick4", "steelblue4"), name="",
                     labels = c("Probability", "Mean Difference")) +
  scale_linetype_manual(values = c(1,2), name="",
                        labels = c("Probability", "Mean Difference")) +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(size = 17, margin=margin(10,0,0,0))) +
  theme(axis.title.y = element_text(size = 17, margin=margin(0,10,0,0), colour = "firebrick4")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15, colour = "firebrick4")) +
  theme(plot.title=element_text(size=18)) +
  theme(plot.margin=unit(c(5.5,18.5,5.5,5.5), "pt")) +
  theme(axis.title.y.right = element_text(vjust=3, colour = "steelblue4")) +
  theme(axis.text.y.right = element_text(colour = "steelblue4"))

## repeat process for case 4
df4 <- read.csv("reg_modes_4.csv")
scale4 = max(df4$V3)

plot4 <- ggplot(df4, aes(x = V1, y = V2)) + theme_bw() +
  geom_line(aes(color = "col1", linetype = "col1"), lwd = 1.25) +
  geom_line(aes(y = V3/scale4, color = "col2", linetype = "col2"), lwd = 1.25) +
  labs(title = bquote("Case 4: ("*italic(u)["1,0"]*","*italic(u)["2,0"]*")"~"= ("*0.5*","*0.99*")")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(sec.axis = sec_axis(~.*scale4, breaks = c(0, 0.04, 0.08, 0.12), name=bquote('|'*italic(D)[2] - italic(D)[1]*'|'))) +
  labs(x=bquote('\n'~italic(n)), 
       y=bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')')) +
  scale_color_manual(values = c("firebrick4", "steelblue4"), name="",
                     labels = c("Probability", "Mean Difference")) +
  scale_linetype_manual(values = c(1,2), name="",
                        labels = c("Probability", "Mean Difference")) +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(size = 17, margin=margin(10,0,0,0))) +
  theme(axis.title.y = element_text(size = 17, margin=margin(0,10,0,0), colour = "firebrick4")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15, colour = "firebrick4")) +
  theme(plot.title=element_text(size=18)) +
  theme(plot.margin=unit(c(5.5,18.5,5.5,5.5), "pt")) +
  theme(axis.title.y.right = element_text(vjust=3, colour = "steelblue4")) +
  theme(axis.text.y.right = element_text(colour = "steelblue4"))

## repeat process for case 5
df5 <- read.csv("reg_modes_5.csv")
scale5 = max(df5$V3)

plot5 <- ggplot(df5, aes(x = V1, y = V2)) + theme_bw() +
  geom_line(aes(color = "col1", linetype = "col1"), lwd = 1.25) +
  geom_line(aes(y = V3/scale5, color = "col2", linetype = "col2"), lwd = 1.25) +
  labs(title = bquote("Case 5: ("*italic(u)["1,0"]*","*italic(u)["2,0"]*")"~"= ("*0.85*","*0.9*")")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(sec.axis = sec_axis(~.*scale5, name=bquote('|'*italic(D)[2] - italic(D)[1]*'|'))) +
  labs(x=bquote('\n'~italic(n)), 
       y=bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')')) +
  scale_color_manual(values = c("firebrick4", "steelblue4"), name="",
                     labels = c("Probability", "Mean Difference")) +
  scale_linetype_manual(values = c(1,2), name="",
                        labels = c("Probability", "Mean Difference")) +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(size = 17, margin=margin(10,0,0,0))) +
  theme(axis.title.y = element_text(size = 17, margin=margin(0,10,0,0), colour = "firebrick4")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15, colour = "firebrick4")) +
  theme(plot.title=element_text(size=18)) +
  theme(plot.margin=unit(c(5.5,18.5,5.5,5.5), "pt")) +
  theme(axis.title.y.right = element_text(vjust=3, colour = "steelblue4")) +
  theme(axis.text.y.right = element_text(colour = "steelblue4"))

## repeat process for case 6; add a legend to this plot
df6 <- read.csv("reg_modes_6.csv")
scale6 = max(df6$V3)

plot6 <- ggplot(df6, aes(x = V1, y = V2)) + theme_bw() +
  geom_line(aes(color = "col1", linetype = "col1"), lwd = 1.25) +
  geom_line(aes(y = V3/scale6, color = "col2", linetype = "col2"), lwd = 1.25) +
  labs(title = bquote("Case 6: ("*theta["1,0"]*","*theta["2,0"]*")"~"= ("*-5*","*8*")")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(sec.axis = sec_axis(~.*scale6, name=bquote('|'*italic(D)[2] - italic(D)[1]*'|'))) +
  labs(x=bquote('\n'~italic(n)), 
       y=bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')')) +
  scale_color_manual(values = c("firebrick4", "steelblue4"), name=NULL,
                     labels = c(bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')'), bquote('|'*italic(D)[2] - italic(D)[1]*'|')),
  ) +
  scale_linetype_manual(values = c(1,2), name=NULL,
                        labels = c(bquote(italic(Pr)*'('*italic(D)[2] <= italic(D)[1]*')'), bquote('|'*italic(D)[2] - italic(D)[1]*'|')),
  ) +
  theme(axis.title.x = element_text(size = 17, margin=margin(10,0,0,0))) +
  theme(axis.title.y = element_text(size = 17, margin=margin(0,10,0,0), colour = "firebrick4")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15, colour = "firebrick4")) +
  theme(plot.title=element_text(size=18)) +
  theme(plot.margin=unit(c(5.5,18.5,5.5,5.5), "pt")) +
  theme(axis.title.y.right = element_text(vjust=3, colour = "steelblue4")) +
  theme(axis.text.y.right = element_text(colour = "steelblue4")) +
  theme(
    legend.position = c(0.97, 0.98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.direction="vertical",
    legend.margin = margin(0,0,0,0)
  ) 

## combine subplots into two columns of three rows each
col1 <- plot_grid(plot1 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.9,0.25,0.25),"cm")),
                  plot3 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.9,0.25,0.25),"cm")),
                  plot5 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.9,0.25,0.25),"cm")),
                  nrow = 3, rel_widths = c(1,1,1))

col2 <- plot_grid(plot2 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.4,0.25,0.75),"cm")),
                  plot4 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.4,0.25,0.75),"cm")),
                  plot6 + theme(plot.margin=unit(c(0.25,0.4,0.25,0.75),"cm")),
                  nrow = 3, rel_widths = c(1,1,1))

fig4 <- plot_grid(col1, col2, ncol = 2, rel_heights = c(1,1))

pdf(file = "Fig4_PR.pdf",   # The directory you want to save the file in
    width = 12.5, # The width of the plot in inches (12.41)
    height = 10) # The height of the plot in inches (10.7)

fig4

dev.off()