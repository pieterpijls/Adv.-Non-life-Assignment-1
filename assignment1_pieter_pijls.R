# Pieter Pijls
# r0387948
# Assignment 1

#load packages
library(tidyverse)
library(ggthemes)
library(plotly)
library(ggplot2)
library(knitr)
library(statmod) 
library(actuar)

#load data
#file.choose()

#Part 1
#------

#Import data
#------------
df <- read.table(file="/Users/pieterpijls/Documents/KU LEUVEN/MFAE/ADVANCED NON LIFE/Assignment 1/SeverityCensoring.txt", header=TRUE,sep=" ")
setwd("/Users/pieterpijls/Documents/KU LEUVEN/MFAE/ADVANCED NON LIFE/Assignment 1/")
x <- df$claimAmount; rc <- df$rc; d <- df$deductible

#Exploratory analysis
#--------------------
summary(df$claimAmount)
summary(df$rc)

#histogram claimAmount
p <- ggplot(df, aes(x = claimAmount), fill = rc) +
  geom_histogram(bins = 10) +
  scale_fill_hc() +
  theme_hc()
ggplotly(p)

#Question 1: Exponantial
#-----------------------

# Create indicator for right-censoring
df$rc <- !is.na(df$rc)
x <- df$claimAmount; rc <- df$rc; d <- df$deductible

#Exponantial
#parset up log likelihood function
loglik.exp <- function(par){
  sum(dexp(x[!rc],par,log=T)) + sum(pexp(x[rc],par,log.p=T,lower.tail=F)) - 
    sum(pexp(d,par,log.p=T,lower.tail=F))
}

# optimize log likelihood function
oo <- optimize(loglik.exp,c(0,1),maximum=T)

#parameter value for exp distr
par.exp <- oo$maximum
AIC.exp <- 2-2*oo$objective


#Lognormal
loglik.lnorm <- function(par){
  sum(dlnorm(x[!rc],par[1],par[2],log=T)) + sum(plnorm(x[rc],par[1],par[2],log.p=T,lower.tail=F)) - 
    sum(plnorm(d,par[1],par[2],log.p=T,lower.tail=F))
}

oo <- optim(c(2*log(mean(x))-log(mean(x^2))/2,sqrt(log(mean(x^2))-2*log(mean(x)))),loglik.lnorm,control=list(fnscale=-1))
par.lnorm <- oo$par
AIC.lnorm <- 4-2*oo$value

#Inverse Gaussian
loglik.invg <- function(par){
  sum(dinvgauss(x[!rc],par[1],par[2],log=T)) +
    sum(pinvgauss(x[rc],par[1],par[2],log.p=T,lower.tail=F)) -
    sum(pinvgauss(d,par[1],par[2],log.p=T,lower.tail=F))
}

oo <- optim(c(mean(x),mean(x)^3/var(x)),loglik.invg,control=list(fnscale=-1))
par.invg <- oo$par
AIC.invg <- 4-2*oo$value

#Burr
loglik.burr <- function(par){
  sum(dburr(x[!rc],shape1=par[1],shape2=par[2],rate=par[3],log=T)) +
    sum(pburr(x[rc],shape1=par[1],shape2=par[2],rate=par[3],log.p=T,lower.tail=F)) -
    sum(pburr(d,shape1=par[1],shape2=par[2],rate=par[3],log.p=T,lower.tail=F))
}

oo <- optim(c(1,1,1),loglik.burr,control=list(fnscale=-1))
par.burr <- oo$par
AIC.burr <- 6-2*oo$value

#Question 2
#----------
source("2014-12-16_ME.R")
nrc <- x; nrc[rc] <- NA
fit.ME <- ME_fit(x, nrc, trunclower = 100, M=5, s=3)
theta <- fit.ME$theta 
shape <- fit.ME$shape 
alpha <- fit.ME$alpha
AIC.ME <- fit.ME$AIC

#Question 3:
#----------
source("2014-12-16_ME.R")
nrc <- x; nrc[rc] <- NA
fit.ME <- ME_fit(x, nrc, trunclower = 100, M=5, s=3)
theta <- fit.ME$theta 
shape <- fit.ME$shape 
alpha <- fit.ME$alpha
AIC.ME <- fit.ME$AIC
AIC.ME

#Question 4
#----------
deds <- d ; loss <- x ; full <- rc
fit <- survfit(Surv(deds, loss, full) ~ 1)
plot(fit, mark.time=F, conf.int=F,lwd=2)

library(ggfortify)
autoplot(fit)
ggplotly()
#adjust axis titles
#try to create interactive graph for other models as well

#Question 5
#----------
surv.exp <- function(y) {pexp(y,par.exp,lower.tail=F)/pexp(100,par.exp,lower.tail=F)}
surv.lnorm <- function(y) {plnorm(y,par.lnorm[1],par.lnorm[2],lower.tail=F)/plnorm(100,par.lnorm[1],par.lnorm[2],lower.tail=F)}
surv.invg <- function(y) {pinvgauss(y,par.invg[1],par.invg[2],lower.tail=F)/pinvgauss(100,par.invg[1],par.invg[2],lower.tail=F)}
surv.burr <- function(y) {pburr(y,par.burr[1],par.burr[2],par.burr[3],lower.tail=F)/pburr(100,par.burr[1],par.burr[2],par.burr[3],lower.tail=F)}
surv.ME <- function(y) {ME_cdf(y, theta, shape, alpha, trunclower = 100, lower.tail = FALSE)}
curve(surv.exp(x),from=0,col=2,add=T)
curve(surv.lnorm(x),from=0,col=3,add=T)
curve(surv.invg(x),from=0,col=4,add=T)
curve(surv.burr(x),from=0,col=5,add=T)
curve(surv.ME(x),from=0,col=6,add=T)
legend('topright', legend = c('KM estimate', 'exp', 'lnorm', 'invgauss', 'burr', 'ME'), col = 1:6, lwd = c(2,1,1,1,1,1))

 curve(surv.burr(x),from=0,col=5,add=T) + curve(surv.ME(x),from=0,col=6,add=T)

library(ggfortify)
autoplot(fit,surv.exp)
ggplotly()

#Question 6
#----------
c(AIC.exp, AIC.lnorm, AIC.invg, AIC.burr, AIC.ME)

#-------
# Part 2
#-------

#Question 7
#----------
loss <- read.table("SecuraRe.txt",header=T)$Loss
n <- length(loss)
k <- 95
sh <- 1200000
thr <- sort(loss)[n-k]
loglik <- function(par) {
  lambda <- par[1]
  alpha <- par[2]
  # indicator
  I <- loss <= thr
  # Likelihood contributions
  L <- I * (n-k)/n * lambda * exp(-lambda*(loss-sh))/(1 - exp(-lambda*(thr-sh))) + 
    (1-I) * k/n * alpha * (loss+1)^(-alpha-1)/(thr+1)^(-alpha)
  # Log Likelihood
  sum(log(L)) 
}
oo <- optim(c(1/mean(loss), log(2)/log(median(loss))), loglik, control=list(fnscale=-1))
lambda <- oo$par[1]
alpha <- oo$par[2]

#Question 8
#----------
plot(ecdf(loss), do.points = FALSE, xlab = 'Claim size', ylab = 'CDF', main = 'Comparison of Empirical CDF and fitted CDF with splicing', xlim = c(sh, max(loss)), lwd = 2)
# Fitted CDF
curve((x >= sh) * ((x <= thr) * (n-k)/n * (1 - exp(-lambda*(x - sh))) / (1 - exp(-lambda*(thr - sh))) + (x > thr) * (1 - k/n * x^(-alpha) / thr^(-alpha))), col=2, lwd=2, add=T)
legend('right', c('Empirical CDF', 'Fitted CDF'), col = c(1, 2), lwd = 2)


