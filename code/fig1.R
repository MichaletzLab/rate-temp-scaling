# Description ----

# Analyses and figures for
# Michaletz, S.T. & Garen, J.C. (in review) On the scaling of biological rates with temperaure.

# Initialize ----
#--Set working directory
setwd("E:/Documents/MS - working/2023 Michaletz Garen temp scaling/temp response data/")

#--Libraries
library(nls.multstart)

#--Load data from Neumeyer (1995)
neu <- read.csv("Data from David Ratkowsky/Neumeyer_1995.csv",header=T)
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
confint(fit_neu1)
# 95% CI
summary(fit_neu1)$parameters[2,1] - qnorm(0.975) * summary(fit_neu1)$parameters[2,2]; summary(fit_neu1)$parameters[2,1] + qnorm(0.975) * summary(fit_neu1)$parameters[2,2]
# Calculate S-S estimates of rate for plotting
predict_SS1 <- data.frame(Temp_K = neu$Temp_K, rate = predict(fit_neu1, newdata=neu$Temp_K))

neu1 <- ggplot () + 
  geom_point(data = neu, aes(x = Temp_K, y = rate), shape = 19, size = 2.75, col = "gray85") +                          
  #geom_line(data = predict_SS0, aes(x = Temp_K, y = rate), color = 'black') + # Plot Sharpe-Schoolfield w/ high-temp deactivation
  geom_line(data = predict_SS1, aes(x = Temp_K, y = rate), color = 'black', linewidth = 0.75) + # Plot Sharpe-Schoolfield w/ low- and high-temp deactivation
  xlab(expression(paste("Temperature (K)"))) + 
  scale_x_continuous(sec.axis = sec_axis(trans = ~ .-273.15 , name = expression(paste('Temperature'~(degree*C))))) +
  ylab(expression(paste("Population growth rate (cells ", cell^{-1}, " ", time^{-1}, ")"))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
neu1