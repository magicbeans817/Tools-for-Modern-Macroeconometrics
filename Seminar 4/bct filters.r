##########################################
####### Hodrick-Prescott filter ##########
##########################################

gdp <- read.csv("C:/Users/Baxa/Documents/Vyuka/Macroeconometrics/Lecture5SpectrumFilters/code/gdpdata.csv", sep=";")
#View(gdp)

gdpCZ <- gdp[,"CZ"]*1e-3
gdpDE <- gdp[,"DE"]*1e-3

#install.packages("mFilter")
library(mFilter)
log_gdpCZ <- log(gdpCZ)
hpcz <- hpfilter(log_gdpCZ,freq=1600,type="lambda")
# residuals(hpcz)
# fitted(hpcz)
plot(hpcz)

##########################################
########## Baxter-King filter ############
##########################################

# bkcz <- bkfilter(log_gdpCZ,pl=6,pu=32)
bkcz <- bkfilter(log_gdpCZ)
summary(bkcz)
plot(bkcz)

### Some plotting ###

par(mfrow=c(1,1))
plot(log_gdpCZ, type="l")
lines(hpcz$trend, col="red")
lines(bkcz$trend, col="blue")

### Spectra of cycles ###

par(mfrow=c(1,2))
hp_spectrum <- spectrum(hpcz$cycle)
bk_spectrum <- spectrum(bkcz$cycle[4:91])
# only part of the series bkcz$cycle is used so that NAs are skipped
# other workaround:
# good <- complete.cases(log_gdpCZ,hpcz$trend,bkcz$trend)
# bkcz_cycle = bkcz$cycle[good]
# this will create a new series only with valid observations
plot(hp_spectrum$freq, hp_spectrum$spec, type = "l")
plot(bk_spectrum$freq, bk_spectrum$spec, type = "l")

# Other filters that are available in library mFilter

# Butterworth filter (band-pass)
# bwfilter(x,freq=NULL,nfix=NULL,drift=FALSE)
# Christiano-Fitzgerald filter (band-pass)
# cffilter(x,pl=NULL,pu=NULL,root=FALSE,drift=FALSE,
#         type=c("asymmetric","symmetric","fixed","baxter-king","trigonometric"#),
#         nfix=NULL,theta=1)
# Trigoniometric filter, pl and pu are lower and upper cut-off frequencies, set in lengths of the period
# trfilter(x,pl=NULL,pu=NULL)


#######################################
########## Hamilton filter ############
#######################################


#install.packages("neverhpfilter")
library(neverhpfilter)
library(xts)
dates <- seq(as.Date("1995-01-01"), length = 94, by = "quarters")
log_gdpCZ_xts <- xts(log_gdpCZ, order.by=dates)

gdp_HAM <- yth_filter(log_gdpCZ_xts, h= 8, p = 4)
par(mfrow=c(1,1))
plot(gdp_HAM["1995/"][,1:2], grid.col = "white", legend.loc = "topleft", main = "Log of Real GDP and trend", panels = 'lines(gdp_HAM["1995/"][,3], type="h", on=NA)')
