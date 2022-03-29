library(dplyr)
library(mFilter)
library(pdfetch)

install.packages("midasr")
library(midasr)



unemp <- read.csv("DP_LIVE_28032022211205262.csv")
unemp <- unemp[,c("TIME", "Value")] # subset only the important stuff

gdp <- read.csv("CLVMNACSCAB1GQFI.csv")
colnames(gdp) <- c("date", "real_gdp")

gdp <- ts(gdp$real_gdp, frequency = 4, start = c(1990, 1))
unemp <- ts(unemp$Value, frequency = 12, start = c(1990, 1))

gdp_x <- window(gdp, end = c(2021, 4))
unemp_y <- window(unemp, end = c(2021, 10))

gdp_ld <- diff(log(gdp_x))*100
unemp_ld <- diff(log(unemp_y)) * 10

gdp_ld %>% length()
unemp_ld %>% length()

# Graph of the data
plot.ts(gdp_ld, xlab = "Time", ylab = "Percentages", col = "blue", ylim = c(-5, 6))
lines(unemp_ld, col = "red")

beta0 <- midasr::midas_r(gdp_ld ~ midasr::mls(gdp_ld, 1, 1) + midasr::mls(unemp_ld, 3:11, 3, nbeta), start = list(unemp_ld = c(1, 1, 1)))
coef(beta0) #prints the estimated coefficients

























