# Example of MIDAS 

# Forecasting GDP growth with the help of monthly employment growth figures.
# Based on Ghysels - Kvedaras - Zemlys: Journal of Statistical Software, 2016
# See https://www.jstatsoft.org/article/view/v072i04.

# Preliminaries
rm(list = ls())
library(midasr)  # loads the midasr package

# Get the data
data("USqgdp")
data("USpayems")

# Convert to ts and subsampling
y <- window(USqgdp, end = c(2011, 2))
x <- window(USpayems, end = c(2011, 7))

# Calculate the log differences
yg <- diff(log(y))*100
xg <- diff(log(x))*100

# It is necessary to equalize the sample size: time series of GDP starts in 1947, employment already in 1939, at the end the sample is equalized to end in 2014Q3. Empty observations are filled with NAs.
nx <- ts(c(NA, xg, NA, NA), start = start(x), frequency = 12)
ny <- ts(c(rep(NA, 33), yg, NA), start = start(x), frequency = 4)

# Graph of the data
plot.ts(nx, xlab = "Time", ylab = "Percentages", col = 4, ylim = c(-5, 6))
lines(ny, col = 2)

# Set the subsample for estimation from 1985Q1 to 2009Q1.
xx <- window(nx, start = c(1985, 1), end = c(2009, 3))
yy <- window(ny, start = c(1985, 1), end = c(2009, 1))

# Estimate the MIDAS models 

# Here Beta polynomial ("nbeta") is used. If Almon lags prefered - change "nbeta" to "nealmon".
# Frequency alignment achieved via mls command. 
# Syntax: mls(x,k,m) ... x = series, k =lags of of the high-frequency variable; m = frequency ratio between high and low frequency.
# In our case, we have low-frequency variable yy, i.e. k = m = 1. The high-frequency variable has 8 lags, starting from lag 3 which implies k = 3:11; m = 3 because there are 3 months within each quarter.
# With restricted MIDAS, it is necessary to set starting values for each variable with restricted coefficients, this is required by underlying nls function. Implicitly, this defines the number of parameters of the constraint function. When not supplied, unrestricted MIDAS is estimated.
#Choice of these parameters, however, not well documented.

beta0 <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbeta), start = list(xx = c(1.7, 1, 5)))
coef(beta0) #prints the estimated coefficients

# point forecast for the 2009Q2 with values of high-frequency variable x provided
xnew <- window(nx, start = c(2009, 4), end = c(2009, 6))
forecast(beta0, newdata = list(xx = xnew), se=TRUE)
forecast(beta0, newdata = list(xx = rep(NA,3)), se=TRUE) #without supplying new data

# Different MIDAS with Note the nbetaMT, which has non-zero weights.
betan <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbetaMT), start = list(xx = c(2, 1, 5, 0)))
coef(betan)
forecast(betan, newdata = list(xx = xnew), se=TRUE)

## Unrestricted MIDAS
um <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3), start = NULL)
coef(um)
forecast(um, newdata = list(xx = xnew), se=TRUE)

#actual value of y in 2009Q2:
print(yg[249,1])

## Compare out-of-sample forecasting performance of all three models

## Split the data into in-sample and out-of-sample
fulldata <- list(xx = window(nx, start = c(1985, 1), end = c(2011, 6)), yy = window(ny, start = c(1985, 1), end = c(2011, 2)))
insample <- 1:length(yy)
outsample <- (1:length(fulldata$yy))[-insample]

## Calculate the individual forecasts of each of the model and their weighted average
avgf <- average_forecast(list(beta0, betan, um), data = fulldata, insample = insample, outsample = outsample)
sqrt(avgf$accuracy$individual$MSE.out.of.sample)
## Apparently, the unrestricted MIDAS has the lowest MSE

avgf$forecast #forecasts of all three models
avgf$avgforecast #average of forecasts

## Simple plot comparing forecasts
plot(avgf$xout, main = "Forecasts of MIDAS regressions vs. actual data", ylim = c(-1,2))
lines(avgf$forecast[,1], col = 1)
lines(avgf$forecast[,2], col = 2)
lines(avgf$forecast[,3], col = 3)
lines(avgf$avgforecast[,1], col= 4)