# Description ----
# Figure 1 plot and analyses from
# Michaletz, S.T. & Garen, J.C. (in review) On the scaling of biological rates with temperature.

# Initialize ----
#--Libraries
library(nls.multstart)
library(ggplot2)
library(strucchange)
library(tidyverse)
library(lmodel2)
library(ggpmisc)
library(scales)
library(gridExtra)
# Uncomment the following lines to install required packages:
# install.packages("nls.multstart")
# install.packages("ggplot2")
# install.packages("strucchange")
# install.packages("tidyverse")
# install.packages("lmodel2")
# install.packages("ggpmisc")
# install.packages("scales")
# install.packages("gridExtra")

#--Load data from Neumeyer (1995)
neu <- read.csv("https://raw.githubusercontent.com/MichaletzLab/rate-temp-scaling/main/data/Neumeyer_1995.csv",header=T)
neu$Temp_C <- neu$Temp_K - 273.15
neu$invT <- 1/(0.00008617*neu$Temp_K)
neu$neg_invT <- -1/(0.00008617*neu$Temp_K)
neu$rate <- neu$Rate_sqrt^2

# Plot Neumeyer (1995) data ----
# Fit S-S with low- and high-temp deactivation (E = 0.97 eV, E_l = 0.87 eV, E_h = 7.83 eV, Tl = 293.58 K, th = 304.75 K)
# Based on Eq. 2, Molnar et al. (2017) but any version would work
fit_neu1 <- nls_multstart(rate ~ B_Tref * exp(-E/8.617E-5*((1/Temp_K)-(1/298.15))) / (1+exp(El/8.617E-5*((1/Tl)-(1/Temp_K)))+exp(Eh/8.617E-5*((1/Th)-(1/Temp_K)))),
                          data = neu,
                          iter = 1000,
                          start_lower = c(B_Tref = 0.005, E = 0, El = 0, Tl = 273.15, Eh = -10, Th = 300),
                          start_upper = c(B_Tref = 0.025, E = 1, El = 5, Tl = 300, Eh = 0, Th = 310),
                          supp_errors = 'Y',
                          na.action = na.omit,
                          lower = c(B_Tref = 0.005, E = 0, El = 0, Tl = 273.15, Eh = -10, Th = 300))
summary(fit_neu1)
# 95% CI
summary(fit_neu1)$parameters[2,1] - qnorm(0.975) * summary(fit_neu1)$parameters[2,2]; summary(fit_neu1)$parameters[2,1] + qnorm(0.975) * summary(fit_neu1)$parameters[2,2]
# Calculate S-S estimates of rate for plotting
predict_SS1 <- data.frame(Temp_K = neu$Temp_K, rate = predict(fit_neu1, newdata=neu$Temp_K))

# Figure 1a: Plot in absolute space
neu1 <- ggplot () + 
  geom_point(data = neu, aes(x = Temp_K, y = rate), shape = 19, size = 2.75, col = "gray85") +                          
  geom_line(data = predict_SS1, aes(x = Temp_K, y = rate), color = 'black', linewidth = 0.75) + # Plot Sharpe-Schoolfield w/ low- and high-temp deactivation
  xlab(expression(paste("Temperature (K)"))) + 
  scale_x_continuous(sec.axis = sec_axis(trans = ~ .-273.15 , name = expression(paste('Temperature'~(degree*C))))) +
  ylab(expression(paste("Population growth rate (cells ", cell^{-1}, " ", time^{-1}, ")"))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
neu1

# Estimate E using piecewise OLS & RMA regression with 2 breakpoints with strucchange() for T < Tpeak: E_1 = 0.86 (0.80 to 0.93), E_2 = 0.44 (0.39 to 0.49)
# See https://stackoverflow.com/questions/70141883/difficulty-fitting-piecewise-linear-data-in-r
# First, truncate 'neu' dataframe to T < Tpeak
neu_trunc <- subset(neu, neu$neg_invT <= neu$neg_invT[which.max(neu$rate)])
# Get a segment size as a fraction of the number of observations
n <- nrow(neu_trunc)
segmts <- 3
h <- (segmts + 1)/n
# Now estimate the breakpoints
b <- breakpoints(log(rate) ~ neg_invT, h = h, breaks = (segmts - 1L), data = neu_trunc)
bp <- neu_trunc[b$breakpoints, "neg_invT"]
bp <- c(bp, Inf)
neu_trunc$grp <- findInterval(neu_trunc$neg_invT, bp, left.open = TRUE)
# Get OLS slopes E and 95% CI from OLS: E_1 = 0.86 (0.80 to 0.93), E_2 = 0.44 (0.39 to 0.49)
neu_trunc %>% group_by(grp) %>% group_map(~ broom::tidy(lm(log(rate) ~ neg_invT, data = .x), conf.int=T))
# Get RMA slopes E and 95% CI from RMA: E_2 = 0.87 (0.81 to 0.94), E_2 = 0.46 (0.42 to 0.52)
neu_trunc %>% group_by(grp) %>% group_map(~ broom::tidy(lmodel2(log(rate) ~ neg_invT, "interval", "interval", data = .x, nperm=99), conf.int=T))

# Fit S-S with low- and high-temp deactivation
# Based on Eq. 2, Molnar et al. (2017) but any version would work
fit_neu2 <- nls_multstart(rate ~ B_Tref * exp(-E/8.617E-5*((1/(-1/(8.617E-5*neg_invT)))-(1/298.15))) / (1+exp(El/8.617E-5*((1/Tl)-(1/(-1/(8.617E-5*neg_invT)))))+exp(Eh/8.617E-5*((1/Th)-(1/(-1/(8.617E-5*neg_invT)))))),
                          data = neu,
                          iter = 1000,
                          start_lower = c(B_Tref = 0.005, E = 0, El = 0, Tl = 273.15, Eh = -10, Th = 300),
                          start_upper = c(B_Tref = 0.025, E = 1, El = 5, Tl = 300, Eh = 0, Th = 310),
                          supp_errors = 'Y',
                          na.action = na.omit,
                          lower = c(B_Tref = 0.005, E = 0, El = 0, Tl = 273.15, Eh = -10, Th = 300))
summary(fit_neu2)
# Calculate S-S estimates of rate for plotting
predict_SS2 <- data.frame(neg_invT = neu$neg_invT, rate = predict(fit_neu2, newdata=neu$neg_invT))

# Estimate Eq (Pawar et al. 2016)
# Truncate 'neu' dataframe to T < Tpeak
neu_trunc2 <- subset(neu, neu$neg_invT <= neu$neg_invT[which.max(neu$rate)])
# Fit Pawar Quadratic model to T < Topt
neu_trunc2 = neu_trunc2 %>% mutate(lnRate = log(rate)) %>% mutate(invT2 = invT*invT)
fit_neu3 = lm(lnRate ~ invT +invT2, data = neu_trunc2)
# Extract Eq
kb <- 8.62e-5
eq = -coef(fit_neu3)[2] - 2*coef(fit_neu3)[3]*(1/(kb*mean(neu_trunc2$Temp_K))) # Eqn. 10, Pawar et al. (2016)
eq
# Calculate Pawar et al. (2016) estimates of rate for plotting
predict_SS3 <- data.frame(neg_invT = neu_trunc2$neg_invT, lnRate = predict(fit_neu3, newdata=data.frame(invT=neu_trunc2$invT, invT2=neu_trunc2$invT2)))
predict_SS3$rate <- exp(predict_SS3$lnRate)

# Figure 1b: Plot all fits in modified Arrhenius space
neu2 <- ggplot() + 
  geom_point(data = neu, aes(x = neg_invT, y = rate), shape = 19, size = 2.75, col = "gray85") +
  geom_line(data = predict_SS2, aes(x = neg_invT, y = rate), color = 'black', linewidth = 0.75) + # Plot Sharpe-Schoolfield with low- and high-temp deactivation
  stat_ma_line(data=subset(neu, neu$neg_invT <= neu$neg_invT[which.max(neu$rate)]), 
               aes(x = neg_invT, y = rate), method="RMA", range.y = "interval", range.x = "interval", color = "#FF00FF", linewidth = 0.75, se=F) + #Plot RMA Arrhenius fit to data below Tpeak
  geom_line(data = predict_SS3, aes(x = neg_invT, y = rate), color = '#FFA500', linewidth = 0.75) + # Plot Pawar et al. (2016) fit
  stat_ma_line(data=neu_trunc, mapping = aes(x = neg_invT, y = rate, group=grp), method="RMA", range.y = "interval", range.x = "interval", color = "#00BFFF", linewidth = 0.75, se=F) + #Plot strucchange piecewise RMA Arrhenius fit to data below Tpeak
  xlab(expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste("Population growth rate (cells ", cell^{-1}, " ", time^{-1}, ")"))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (-1/(.*0.00008617))-273.15 , name = expression(paste('Temperature'~(degree*C))))) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
neu2

# Figure 1: Plot panels together ----
grid.arrange(neu1, neu2, ncol=2)



