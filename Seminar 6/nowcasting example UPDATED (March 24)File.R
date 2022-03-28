# Nowcasting using factor model and structural time series model

########## Factor model - Illustration ################

# Factor model suggested by Giannone et al. 2008 as implemented in nowcasting package.
# https://github.com/nmecsys/nowcasting
# This particular application is based on NY Fed Nowcasting report
# Model description at https://www.newyorkfed.org/research/staff_reports/sr830
# The goal is to estimate the GDP growth in 2016Q4 and make forecast for 2017, with the knowledge of flash GDP estimate for 2016Q4 and some surveys from early January 2017.

rm(list = ls())
library(nowcasting)

data(NYFED)

# organizing the dataset
base <- NYFED$base #contains the main data
View(base) #note the ragged edge at the last observations. The last observation is for January 2017:
NYFED$Time

blocks <- NYFED$blocks$blocks 
# assigns whether the series belongs to one of four categories
trans <- NYFED$legend$Transformation 
# defines transformation to achieve stationarity for all variables
frequency <- NYFED$legend$Frequency #defines frequency of all variables

# Bpanel transforms original monthly time series to stationary representations, NAs and outliers can be replaced.
x <- Bpanel(base = base, trans = trans, NA.replace = F, na.prop = 1)

# Nowcasting: projection from common factors
nowEM <- nowcast(formula = GDPC1 ~ ., data = x, r = 1, p = 1, 
                 method = "EM", blocks = blocks, frequency = frequency) 
# selected estimation method, the number r of dynamic factors (for each block of variables), the lag order of the factors p

# The function nowcast returns object yfcst - block of original series, in-sample estimates, out-of-sample forecasts + estimated dynamic factors + forecasts of exogenous variables.
# In-sample predictions till 2016Q4 (note that for the last quarter not all variables have observation); out-of-sample predictions for 2017Q1 - 2017Q4.

nowEM$yfcst #printing original series, in-sample forecasts, out-of-sample forecasts. 

nowcast.plot(nowEM)
#nowcast.plot(nowEM, type = "fcst") #equivalent to nowcast.plot(nowEM)
nowcast.plot(nowEM, type = "factors")


########## Basic structural time series model ################

# Alternative - bsts package - bayesian estimation of structural time series model
# https://www.unofficialgoogledatascience.com/2017/07/fitting-bayesian-structural-time-series.html
# This model fits the structural time series model of short-term forecasts of weekly claims of unemployment insurance in the U.S. with the help of google trends data that are available sooner. It is shown that inclusion of google trends data improves forecasting performance.

library(bsts)     # load the bsts package
data(iclaims)     # bring the initial.claims data into scope

# Building the structural time series model from its components
ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)

# Estimate the structural time series model of series iclaimsNSA, 1000 iterations of MCMC algorithm are used.

model1 <- bsts(initial.claims$iclaimsNSA,
               state.specification = ss,
               niter = 1000)

plot(model1)
plot(model1, "components")

pred1 <- predict(model1, horizon = 12)
plot(pred1, plot.original = 156)

# Now, Google trends data are added.
# First, fit a bsts model with expected model size 1, the default.
# The algorithm sets the initial probability to each coefficients zero, and updates the probabilities to get posterior. Variable with largest posterior density is included in the final model.
model2 <- bsts(iclaimsNSA ~ .,
               state.specification = ss,
               niter = 1000,
               data = initial.claims)

plot(model2)
plot(model2, "components")
plot(model2, "coef")

# Fit a bsts model with expected model size 5, to include more coefficients.
model3 <- bsts(iclaimsNSA ~ .,
               state.specification = ss,
               niter = 1000,
               data = initial.claims,
               expected.model.size = 5)

# Model diagnostics: Did the Google data help?
# This function compares one-step ahead prediction errors.

CompareBstsModels(list("Model 1" = model1,
                       "Model 2" = model2,
                       "Model 3" = model3),
                  colors = c("black", "red", "blue"))

# The bottom panel shows the original series. The top panel shows the cumulative total of the mean absolute one step prediction errors for each model. The final time point in the top plot is proportional to the mean absolute prediction error for each model, but plotting the errors as a cumulative total lets you see particular spots where each model encountered trouble, rather than just giving a single number describing each model???s predictive accuracy. Figure  shows that the Google data help explain the large spike near 2009, where model 1 accumulates errors at an accelerated rate, but models 2 and 3 continue accumulating errors at about the same rate they had been before. The fact that the lines for models 2 and 3 overlap means that the additional predictors allowed by the relaxed prior used to fit model 3 do not yield additional predictive accuracy.