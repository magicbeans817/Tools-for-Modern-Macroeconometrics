# Macroeconometrics
# Class 5: State space models
# Jaromir Baxa, IES FSV UK

# Local level model and local trend model with StructTS and DLM
# Source: https://kevinkotze.github.io/ts-4-tut/ and Petris, G., & Petrone, S. (2011). State Space Models in R. Journal of Statistical Software, 41(4), 1–25. https://doi.org/10.18637/jss.v041.i04 (https://www.jstatsoft.org/article/view/v041i04)

rm(list = ls())
graphics.off()
#devtools::install_github("KevinKotze/tsm")
#install.packages("dlm")
library(tsm)
library(dlm)
library(zoo)
library(forecast)

# Get the data from TSM package, derive GDP deflator from nominal and real GDP
dat <- sarb_quarter$KBP6006L/sarb_quarter$KBP6006D
dat.tmp <- diff(log(na.omit(dat)) * 100, lag = 1)
head(dat)

# create time series object
inf <- ts(dat.tmp, start = c(1960, 2), frequency = 4)
plot.ts(inf)

######StructTS (base)#########

# StructTS function (base package)
# StructTS(x, type = c("level", "trend", "BSM"), init = NULL, fixed = NULL, optim.control = NULL)
# BSM - includes seasonal terms based on frequency of the data, fixed is for coefficients that shall be constant.
model1 <- StructTS(inf, type="level")
model1
model2 <- StructTS(inf, type="trend")
model2

# plotting filtered and smoothed series
plot(inf, type = "o", main="Inflation: Local level (red) and local trend (blue) models")
lines(fitted(model1), lty = "dotted", lwd = 3, col = "red") 
lines(tsSmooth(model1), lwd = 4, col = "red")
# model 2 has two elements - level and slope
model2sm <- tsSmooth(model2)
lines(fitted(model2)[,1], lty = "dotted", lwd = 3, col = "blue")
lines(model2sm[,1], lwd = 4, col = "blue")

# showing decomposition to trend, slope and seasonals
model3 <- StructTS(inf, type = "BSM")
model3_decomp <- cbind(inf,fitted(model3))
colnames(model3_decomp) <- c("data","level","slope", "seasonal")
plot(model3_decomp, main="Decomposition of Inflation")

# forecasting
par(mfrow=c(2,1))
model1_f <- forecast(model1)
plot(model1_f)
model2_f <- forecast(model2)
plot(model2_f)

######DLM package#########

# The StructTS fairly reliable and straightforward, but the options are somewhat restrictive.
# More complex models, i.e. with additional regressors, multivariate models or estimation with Bayesian techniquest can be done in DLM package.
# Bayesian estimation of simple structural models can be, however, estimated using the bsts package. https://rdrr.io/cran/bsts/man/

# Local level model using dlmModPoly; order = 1 => local level, order = 2 => local trend model
# Next, we set parameters for the errors in measurement and state equation
dlm1 <- function(parm) {
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}

# Fitting the model: parm - setting initial values, here zeros; build - uses previous object. 
dlm1_fit <- dlmMLE(inf, parm = c(0,0), build = dlm1, hessian = TRUE)
(conv <- dlm1_fit$convergence)

# Get some statistics
loglik <- dlmLL(inf, dlmModPoly(1))
n.coef <- 2
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(inf))) * (n.coef)
dlm1_mod <- dlm1(dlm1_fit$par)
obs.error.var <- V(dlm1_mod)
state.error.var <- W(dlm1_mod)

# Filtering and smoothing
dlm1_filt <- dlmFilter(inf,dlm1_mod)
dlm1_smooth <- dlmSmooth(dlm1_filt)

# Get residuals and create time series object from smoothed series 
resids <- residuals(dlm1_filt, sd = FALSE)
mu <- dropFirst(dlm1_smooth$s)
mu.1 <- mu[1]
mu.end <- mu[length(mu)]

# Diagnostics
shapiro.test(resids)
par(mfrow = c(1,1))
hist(resids, prob = TRUE, breaks = seq(-4.5,7, length.out = 30))
par(mfrow = c(1,2))
acf(resids)
pacf(resids)
Box.test(resids,lag=12,type = "Ljung")

# Confidence intervals

conf.tmp <- unlist(dlmSvd2var(dlm1_smooth$U.S, dlm1_smooth$D.S))
conf <- ts(as.numeric(conf.tmp)[-1], start = c(1960, 2), 
           frequency = 4)
wid <- qnorm(0.05, lower = FALSE) * sqrt(conf)

conf.pos <- mu + wid
conf.neg <- mu - wid

par(mfrow = c(1, 1))
plot.ts(inf, col = "black", xlab = "", ylab = "", lwd = 1)
lines(mu, col = "red", lwd = 3)
lines(conf.pos, col = "darkgrey", lty = "dashed", lwd = 2)
lines(conf.neg, col = "darkgrey", lty = "dashed", lwd = 2)
legend("topright", legend = c("Observed Deflator", "Stochastic level", "Confidence Interval"), lwd = c(1, 2, 2), col = c("black", "red", "darkgrey"), bty = "n")

# One can estimate Local trend model in a similar way. dlmModPoly needs to be of order 2 and one more parameter needs to be estimated.

dlm2 <- function(parm) {
  dlmModPoly(order = 2, dV = exp(parm[1]), dW = exp(parm[2:3]))
}

dlm2_fit <- dlmMLE(inf, parm = c(0,0,0), build = dlm2, hessian = TRUE)
(conv <- dlm2_fit$convergence)

dlm2_mod <- dlm2(dlm2_fit$par)
dlm2_filtered <- dlmFilter(inf, dlm2_mod)
dlm2_smoothed <- dlmSmooth(dlm2_filtered)

resids <- residuals(dlm2_filtered, sd = FALSE)
mu <- dropFirst(dlm2_smoothed$s[, 1])
upsilon <- dropFirst(dlm2_smoothed$s[, 2])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
ups.1 <- upsilon[1]
ups.end <- upsilon[length(mu)]

#Plotting the decomposition
par(mfrow = c(3, 1))
plot.ts(inf, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(mu, col = "black")
legend("topright", legend = c("Observed Deflator", "Stochastic level"), lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(upsilon, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
legend("topright", legend = "Slope", lwd = 2, col = "darkgrey", bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", lwd = 2)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 2, col = "darkgrey", bty = "n")

# Plot with confidence intervals

conf.tmp <- unlist(dlmSvd2var(dlm2_smooth$U.S, dlm2_smooth$D.S))
conf <- ts(as.numeric(conf.tmp)[-1], start = c(1960, 2), frequency = 4)
wid <- qnorm(0.05, lower = FALSE) * sqrt(conf)

conf.pos <- mu + wid
conf.neg <- mu - wid

par(mfrow = c(1, 1))
plot.ts(inf, col = "black", xlab = "", ylab = "", lwd = 1)
lines(mu, col = "red", lwd = 3)
lines(conf.pos, col = "darkgrey", lty = "dashed", lwd = 2)
lines(conf.neg, col = "darkgrey", lty = "dashed", lwd = 2)
legend("topright", legend = c("Observed Deflator", "Stochastic level", "Confidence Interval"), lwd = c(1, 2, 2), col = c("black", "red", "darkgrey"), bty = "n")

# Forecast
dlm2_f <- dlmForecast(dlm2_filtered, nAhead = 12)
var.2 <- unlist(dlm2_f$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- cbind(dlm2_f$f, dlm2_f$f + wid.2, dlm2_f$f - wid.2)
comb.state <- cbind(mu, conf.pos, conf.neg)
result <- ts(rbind(comb.state, comb.fore), start = c(1960,2), frequency = 4)

par(mfrow = c(1, 1))
plot.ts(result, col = c("red", "darkgrey", "darkgrey"), plot.type = "single", 
        xlab = "", ylab = "", lty = c(3, 2, 2), ylim = c(-2, 8))
lines(inf, col = "black", lwd = 1.5)
abline(v = c(2017, 1), col = "blue", lwd = 1, lty = 3)
legend("topleft", legend = c("Observed Deflator", "Stochastic level"), 
       lwd = c(1.5, 1), col = c("darkgrey", "black"), bty = "n")

######### Missing values ######

# Appealing property of state space models is that they can handle series with missing values quite easily: they are treated as if they were forecasts.

inf.mis <- inf
inf.mis[70:82] <- NA
plot(inf.mis)

dlm2mis <- function(parm) {
  dlmModPoly(order = 2, dV = exp(parm[1]), dW = exp(parm[2:3]))
}

dlm2mis_fit <- dlmMLE(inf.mis, parm = c(0,0,0), build = dlm2mis, hessian = TRUE)
(conv <- dlm2mis_fit$convergence)

dlm2mis_mod <- dlm2mis(dlm2mis_fit$par)
dlm2mis_filtered <- dlmFilter(inf.mis, dlm2mis_mod)
dlm2mis_smoothed <- dlmSmooth(dlm2mis_filtered)

resids <- residuals(dlm2mis_filtered, sd = FALSE)
mu <- dropFirst(dlm2mis_smoothed$s[, 1])
upsilon <- dropFirst(dlm2mis_smoothed$s[, 2])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
ups.1 <- upsilon[1]
ups.end <- upsilon[length(mu)]

# Plotting the decomposition
par(mfrow = c(3, 1))
plot.ts(inf.mis, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
lines(mu, col = "black")
legend("topright", legend = c("Observed Deflator", "Stochastic level"), lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(upsilon, col = "darkgrey", xlab = "", ylab = "", lwd = 2)
legend("topright", legend = "Slope", lwd = 2, col = "darkgrey", bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", lwd = 2)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 2, col = "darkgrey", bty = "n")

# Plot with confidence intervals

conf.tmp <- unlist(dlmSvd2var(dlm2mis_smooth$U.S, dlm2mis_smooth$D.S))
conf <- ts(as.numeric(conf.tmp)[-1], start = c(1960, 2), frequency = 4)
wid <- qnorm(0.05, lower = FALSE) * sqrt(conf)

conf.pos <- mu + wid
conf.neg <- mu - wid

par(mfrow = c(1, 1))
plot.ts(inf.mis, col = "black", xlab = "", ylab = "", lwd = 1)
lines(mu, col = "red", lwd = 3)
lines(conf.pos, col = "darkgrey", lty = "dashed", lwd = 2)
lines(conf.neg, col = "darkgrey", lty = "dashed", lwd = 2)
legend("topright", legend = c("Observed Deflator", "Stochastic level", "Confidence Interval"), lwd = c(1, 2, 2), col = c("black", "red", "darkgrey"), bty = "n")
