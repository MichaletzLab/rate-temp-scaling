# Description ----
# Figure 2 plot and analyses from
# Michaletz, S.T. & Garen, J.C. (in review) On the scaling of biological rates with temperature.

# Code to reproduce photosynthesis model and activation energy
# calculation from Allen et al. (2005). Main functions are:
#   v_chlo - gives photosynthetic rate as a function of temperature (K)
#            using the parameters given in Bernacchi et al. (2001)
#   act_e  - can be used to reproduce the activaiton energy value given
#            in Allen
#
#            Josef Garen - 14 November 2018
# Modified by Sean Michaletz (sean.michaletz@gmail.com), 14 September 2023

# Initialize ----
#--Libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggtext)
library(nls.multstart)
# Uncomment the following lines to install required packages:
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("ggtext")
# install.packages("nls.multstart")

# Maximum rate of rubisco carboxylation per chloroplast
# function of temperature (K)
# return value in micromol m^-2 s^-1
V_c_max <- function(temp_K) {
  c   <- 26.35           # dimensionless scaling constant
  E_a <- 65.33           # activation energy (KJ/mol)
  R   <- 8.31446*(10^-3) # molar gas constant (KJ/mol*K)
  
  return(exp(c - (E_a/(R*temp_K))))
}

# Maximum rate of rubisco oxygenation (photorespiration) per chloroplast
# function of temperature (K)
# return value in micromol m^-2 s^-1
# Note Bernacchi claims this function is normalized at T=25C, but not true
V_o_max <- function(temp_K) {
  c   <- 22.98           # dimensionless scaling constant
  E_a <- 60.11           # activation energy (KJ/mol)
  R   <- 8.31446*(10^-3) # molar gas constant (KJ/mol*K)
  
  return(exp(c - (E_a/(R*temp_K))))
}

# Michaelis-Menton constant for CO2
# function of temperature (K)
# Return value in micromol/mol (note: different from K_o)
K_c <- function(temp_K) {
  c   <- 38.05           # dimensionless scaling constant
  E_a <- 79.43           # activation energy (KJ/mol)
  R   <- 8.31446*(10^-3) # molar gas constant (KJ/mol*K)
  
  return(exp(c - (E_a/(R*temp_K))))
}

# Michaelis-Menton constant for O2
# function of temperature (K)
# Return value in millimol/mol (note: different from K_o)
K_o <- function(temp_K) {
  c   <- 20.3            # dimensionless scaling constant
  E_a <- 36.38           # activation energy (KJ/mol)
  R   <- 8.31446*(10^-3) # molar gas constant
  
  return(exp(c - (E_a/(R*temp_K))))
}

# Ratio of rubisco oxidation to carboxylation
# function of temperature (K)
# Returns unitless ratio
phi <- function(temp_K) {
  pres_o <- 210  # O2 concentration (millimol/mol)
  pres_c <- 280  # CO2 concentration (micromol/mol)
  
  num = V_o_max(temp_K)*K_c(temp_K)*pres_o   # compute numerator
  denom = V_c_max(temp_K)*K_o(temp_K)*pres_c # compute denominator
  
  return(num/denom)
}

# Rate of photosynthesis per chloroplast
# funciton of temperature (K)
# return value in micromol m^-2 s^-1
v_chlo <- function(temp_K) {
  pres_o <- 210 # O2 concentration (millimol/mol)
  pres_c <- 280 # CO2 concentration (millimol/mol)
  
  num = (1-(phi(temp_K)/2))*V_c_max(temp_K)                 # compute numerator
  denom = 1+((K_c(temp_K)/pres_c)*(1+(pres_o/K_o(temp_K)))) # compute denominator
  
  return(num/denom)
}

# Compute activation energy according to the method given in Allen
# Accepts high and low temperatures (K) as args
# Finds ratio of photosynthesis at each temp and fits a boltzmann curve using
# only the extreme values
# return value in eV
act_e <- function (T_low, T_high) {
  k = 8.617e-5 # Boltzmann constant in eV K^-1
  
  return(k*log(v_chlo(T_high)/v_chlo(T_low))/((1/T_low)-(1/T_high)))
}

const_c <- function (T_low, T_high) {
  k = 8.617e-5 # Boltzmann constant in eV K^-1
  
  Ea = k*log(v_chlo(T_high)/v_chlo(T_low))/((1/T_low)-(1/T_high))
  
  return(log(v_chlo(T_low)*exp(Ea/(k*T_low))))
  
}

# Generate a plot which illustrates Allen's method of estimating Ea of photosynthesis
k = 8.617e-5
temp_K = 273:323
sim_photo = v_chlo(temp_K)


Ea = act_e(273,303)
cc = const_c(273,303)
BA_photo = exp(cc - Ea/(k*temp_K))

df = data.frame(temperature = temp_K-273.15,
                Photo_FvCB = sim_photo, 
                Photo_BAfit = BA_photo)

df = df %>% mutate(Photo_BAfit = Photo_BAfit/max(Photo_FvCB),
                   Photo_FvCB = Photo_FvCB/max(Photo_FvCB))

# Figure 2: Plot A-T in modified Arrhenius space----
df$neg_invT <- -1/(0.00008617*(df$temperature + 273.15))
fig2 <- ggplot() + 
  geom_point(data = df, aes(x = neg_invT, y = Photo_FvCB), shape = 19, size = 1, col = "#00B050") +
  geom_line(data = df, aes(x = neg_invT, y = Photo_BAfit), color = 'black', linewidth = 0.75, lty = 2) +
  geom_point(data = df[c(1,31),], aes(x = neg_invT, y = Photo_FvCB), color = "#0000FF",fill = "#0000FF", pch = 21, size = 2.75) +
  xlab(expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste("Photosynthesis (relative units)"))) +
  annotate(geom = "richtext", label = "<span style='color: #00B050;'>FvCB prediction</span><br><span style='color: #0000FF;'>Fitted data </span><br>Arrhenius fit",
           x = -1/(0.00008617*(32+273.15)), y = 0.375, fill = NA, label.color = NA) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (-1/(.*0.00008617))-273.15 , name = expression(paste('Temperature'~(degree*C))))) +
  scale_y_continuous(trans="log", labels = scales::number_format(accuracy = 0.01), limits = c(0.25,1.25)) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
fig2

# Estimate activation energy and CI from SS fit
michaletz_2021 <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25) { 
  # temp   : Temperature values to evaluate at (C)
  # J_ref  : Rate at T_ref (units depend on rate)
  # E      : Activation energy (eV; E > 0)
  # E_D    : Deactivation energy (eV; E_D > 0)
  # T_opt  : Optimum or peak T (C)
  # T_ref  : Reference temperature for normalization (C)
  
  k <- 8.62e-5            # Boltzmann's constant (eV K-1)
  temp = temp + 273.15    # Convert to Kelvin
  T_ref = T_ref + 273.15  # Convert to Kelvin
  T_opt = T_opt + 273.15  # Convert to Kelvin
  
  # Evaluate and return
  return( J_ref * exp( E * (1/(k*T_ref) - 1/(k*temp)) ) / ( 1 + E/(E_D - E) * exp( (E_D/k) * (1/T_opt - 1/temp) ) ) )
}

df$T_ref = 10
fit <- nls_multstart(Photo_FvCB ~ michaletz_2021(temperature, J_ref, E, E_D, T_opt, T_ref),
                     data = df,
                     iter = 1000,
                     start_lower = c(J_ref = 0, E = 0, E_D = 0, T_opt = 0),
                     start_upper = c(J_ref = 20, E = 2, E_D = 5, T_opt = 50),
                     supp_errors = 'Y',
                     na.action = na.omit,
                     lower = c(J_ref = -10, E = 0, E_D = 0, T_opt = 0))


confint(fit)
