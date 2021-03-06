###########################################
########## Spectra of series ##############
###########################################


# Function that calculates the fast fourier transformation
# and plots the periodogram
fftMy <- function(y,T,t,smoothFact){
  
  # periodogram    ... y in term of scales, x in term of periods
  # vs. spectogram ... commented lines 
  
  Ft <- abs(fft(y)/(T/2))
  #Ft <- abs(Ft)^2 # y axis scaled
  Ft <- Ft[1:ceiling(T/2)]
  
  #freqs <- t - 2*pi/T
  #freqs <- freqs[1:ceiling(T/2)] # x axis in angular, y axis unscaled
  freqs <- 0:(ceiling(T/2)-1) # x axis in periods, y axis scaled
  
  plot(freqs, Ft, type="h", xlab="Period", ylab=expression(f(x))) #, xlim=c(0, pi))
  points(freqs, Ft, pch=20)
  
  # smoothing -> spectrum
  # smoothFact = T/5 
  lines(smooth.spline(freqs,Ft, df=smoothFact), col="OrangeRed")
  
}

# White noise
yy <- rnorm(100)
par(mfrow=c(1,2))
plot(yy, type = "l")

T <- length(yy)
t <- (1:T)*2*pi/T
fftMy(yy,T,t,50)

# unit root: just add cumsum

zz <- cumsum(rnorm(100))
par(mfrow=c(1,2))
plot(zz, type = "l")

T <- length(zz)
t <- (1:T)*2*pi/T
fftMy(zz,T,t,50)

# MA(1) 
theta <- +0.9
yyE <- rnorm(100)
yy <- yyE[2:100]+theta*yyE[1:100-1]
par(mfrow=c(1,2))
plot(yy, type = "l")
T <- length(yy)
t <- (1:T)*2*pi/T
fftMy(yy,T,t,50)


# AR(1)
sigma <- +0.9
yyE <- rnorm(100)
yy <- yyE
for (jj in 2:100){
  yy[jj] = yy[jj-1]*sigma + yyE[jj]
}
par(mfrow=c(1,2))
plot(yy, type = "l")
T <- length(yy)
t <- (1:T)*2*pi/T
fftMy(yy,T,t,50)

# remember, that MA(1) <=> AR(inf)

# try opposite signs

##########################################
########## Real world example ############
##########################################

gdp <- read.csv("gdpdata.csv", sep=";")
View(gdp)

gdpCZ <- gdp[,"CZ"]*1e-3
gdpDE <- gdp[,"DE"]*1e-3
plot(gdpCZ, type = "l")
lines(gdpDE,col='red')

par(mfrow=c(1,2))

# CZ
plot(gdpCZ, type = "l")
T <- length(gdpCZ)
t <- (1:T)*2*pi/T
fftMy(gdpCZ,T,t,T/5)
#now first differences
d.gdpCZ = diff(log(gdpCZ))*100
plot(d.gdpCZ, type = "l")
T <- length(d.gdpCZ)
t <- (1:T)*2*pi/T
fftMy(d.gdpCZ,T,t,T/5)

# DE
plot(gdpDE, type = "l")
T <- length(gdpDE)
t <- (1:T)*2*pi/T
fftMy(gdpDE,T,t,T/5)
d.gdpDE = diff(log(gdpDE))*100
plot(d.gdpDE, type = "l")
T <- length(d.gdpDE)
t <- (1:T)*2*pi/T
fftMy(d.gdpDE,T,t,T/5)
#now second differences :-)
dd.gdpDE = diff(d.gdpDE)
plot(dd.gdpDE, type = "l")
T <- length(dd.gdpDE)
t <- (1:T)*2*pi/T
fftMy(dd.gdpDE,T,t,T/5)

# Built-in function "spectrum"

par(mfrow=c(1,1))
GDP.spectrum = spectrum(gdpCZ)
# However, the y-axis has a log-scale
par(mfrow=c(1,2))
plot(GDP.spectrum$spec, type="l")
plot(GDP.spectrum$freq, GDP.spectrum$spec, type = "l")
