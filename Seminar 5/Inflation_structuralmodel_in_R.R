library(pdfetch)

# pdfetch_EUROSTAT("namq_gdp_c", FREQ="Q", S_ADJ="NSA", UNIT="MIO_EUR", 
#                 INDIC_NA="B1GM", GEO=c("DE","UK"))

cpi <- pdfetch_EUROSTAT("prc_hicp_midx", FREQ="M", UNIT="I15", COICOP = "CP00", GEO="CZ")

#Converting xts to ts object using library tsbox
library(tsbox)
x <- ts_first_of_period(cpi) # tsbox requires first day of a period in time stamp
cpi_cz <- ts_ts(x)

inflation_mom <- diff(log(cpi_cz) * 100, lag = 1)
plot.ts(inflation_mom)

model1 <- StructTS(inflation_mom, type="level")
model1
model2 <- StructTS(inflation_mom, type="trend")
model2

# plotting filtered and smoothed series
plot(inflation_mom, type = "o", main="Inflation: Local level (red) and local trend (blue) models")
lines(fitted(model1), lty = "dotted", lwd = 3, col = "red") 
lines(tsSmooth(model1), lwd = 4, col = "red")
# model 2 has two elements - level and slope
model2sm <- tsSmooth(model2)
lines(fitted(model2)[,1], lty = "dotted", lwd = 3, col = "blue")
lines(model2sm[,1], lwd = 4, col = "blue")

# showing decomposition to trend, slope and seasonals
model3 <- StructTS(inflation_mom, type = "BSM")
model3_decomp <- cbind(inflation_mom,fitted(model3))
colnames(model3_decomp) <- c("data","level","slope", "seasonal")
plot(model3_decomp, main="Decomposition of Inflation")

# forecasting
model1_f <- forecast(model1)
plot(model1_f)
model2_f <- forecast(model2)
plot(model2_f)
model3_f <- forecast(model3)
plot(model3_f)

#Inflation y-o-y - but here, the MLE converges to variance of measurement equation = 0.

inflation_yoy <- diff(log(cpi_cz) * 100, lag = 12)
plot.ts(inflation_yoy)

model4 <- StructTS(inflation_yoy, type="level", init = c(1,1))
model4
model5 <- StructTS(inflation_yoy, type="BSM")
model5
model5_decomp <- cbind(inflation_yoy,fitted(model5))
colnames(model5_decomp) <- c("data","level","slope", "seasonal")
plot(model5_decomp, main="Decomposition of Inflation")

# plotting filtered and smoothed series
plot(inflation_yoy, type = "o", main="Inflation: Local level model (red) and BSM (blue)")
lines(fitted(model4), lty = "dotted", lwd = 3, col = "red") 
lines(tsSmooth(model4), lwd = 4, col = "red")
model5sm <- tsSmooth(model5)
lines(fitted(model5)[,1], lty = "dotted", lwd = 3, col = "blue")
lines(model5sm[,1], lwd = 4, col = "blue")

# forecasting
model4_f <- forecast(model4)
plot(model4_f)
model5_f <- forecast(model5)
plot(model5_f)