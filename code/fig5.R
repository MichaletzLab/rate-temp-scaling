# Code for figure 5 showing dependence of Ea estimates on noise level
# and temperature range for SS, Arrhenius, and Eq
# Michaletz ST & Garen JC. 2024. Hotter is not (always) better: Embracing unimodal scaling of biological rates with temperature.

library(rTPC)
library(tidyverse)
library(ggplot2)
library(nls.multstart)
library(gridExtra)

# Returns Arrhenius function for the given parameters
arrhenius <- function(temp, J_ref, E, T_ref = 25) { 
  # temp   : Temperature values to evaluate function at (C)
  # J_ref  : Rate at T_ref (units depend on rate)
  # E      : Activation energy (eV; E > 0)
  # T_ref  : Reference temperature for normalization (C)
  
  k <- 8.62e-5           # Boltzmann's constant (eV K-1)
  temp = temp + 273.15   # Convert to Kelvin
  T_ref = T_ref + 273.15 # Convert to Kelvin
  
  # Evaluate and return
  return( J_ref * exp( E * (1/(k*T_ref) - 1/(k*temp)) ) ) 
}

# Returns Sharpe-Schoolfield function for the given parameters
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


# Generate a full TPC
tpc = data.frame(temperature = 0:50,
                 performance = pawar_2018(0:50, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10))

# Compute 1st and 2nd derivatives of tpc
tpc$dpdt = (c(tpc$performance,NA)-c(NA, tpc$performance))[-1]
tpc$d2pdt2 = (c(tpc$dpdt,NA)-c(NA, tpc$dpdt))[-51]


# Loop from 5 degrees to 45 degrees
# Fit Arrhenius over data
k=1
results_all = data.frame()
results_compiled = data.frame()
set.seed(0)

# Set noise level l; loop over 3 levels
for(l in c(0.005,0.01,0.015)) {
  
  # At each noise level, do 50 reps
  for(j in 1:50) {  
    results = data.frame()
    
    # Loop over temperature ranges from 5-45
    for (i in seq(5,45,2)) {
      
      # Generate TPC
      k=k+1
      cur = data.frame(temperature = i*(0:9)/9,
                       performance = pawar_2018( i*(0:9)/9, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10))
      
      # Add noise to TPC
      cur$performance = cur$performance + rnorm(n = 10, mean = 0, sd = l)
      
      # Set reference temperature
      cur$T_ref=10
      
      # Fit Arrhenius
      fit = nls(performance ~ arrhenius(temperature, J_ref, E, T_ref = 10),
                data = cur,
                start = list(J_ref = 1, E = 0.4))
      
      # Fit Sharpe-Schoolfield
      fit2 <- nls_multstart(performance ~ michaletz_2021(temperature, J_ref, E, E_D, T_opt, T_ref),
                            data = cur,
                            iter = 1000,
                            start_lower = c(J_ref = 0, E = 0, E_D = 0, T_opt = 0),
                            start_upper = c(J_ref = 20, E = 2, E_D = 5, T_opt = 50),
                            supp_errors = 'Y',
                            na.action = na.omit,
                            lower = c(J_ref = -10, E = 0, E_D = 0, T_opt = 0))
      
      # Fit Pawar Quadratic model
      kb <- 8.62e-5 
      cur2 = cur %>% mutate(invTemp = 1/(kb*(temperature+273.15)), lnPerf = log(performance)) %>% 
        mutate(invTemp2 = invTemp*invTemp)
      fit3 = lm(lnPerf ~ invTemp +invTemp2, 
                data = cur2)
      # Extract Eq
      eq = -coef(fit3)[2] - 2*coef(fit3)[3]*(1/(kb*mean(cur$temperature+273.15)))

      # Compile Ea estimates
      cur_results = data.frame(Tmin = 0,
                               Tmax = i,
                               E_sch = coef(fit2)[2],
                               E_arr = coef(fit)[2],
                               E_q = eq)
      results = bind_rows(results, cur_results)
      
    }
    results = results %>% mutate(run = j)
    results_all = bind_rows(results_all, results)
    
  }
    
  results_all = results_all %>% mutate(error_level = l)
  results_compiled = bind_rows(results_compiled, results_all)
}

# Compute means and SD
mean_E = results_compiled %>% 
  group_by(Tmax) %>% 
  summarize(mean_E_arr = mean(E_arr), mean_E_sch = mean(E_sch), mean_E_q = mean(E_q)) %>% 
  gather(key = "E_type", value = "E", mean_E_arr, mean_E_sch, mean_E_q)
mean_sd_arr = results_compiled %>% 
  group_by(Tmax, error_level) %>% 
  summarize(sd_E_arr = sd(E_arr)) %>% 
  spread(key = "error_level", value = "sd_E_arr") %>% 
  mutate(E_type = "mean_E_arr")
mean_sd_sch = results_compiled %>% 
  group_by(Tmax, error_level) %>% 
  summarize(sd_E_sch = sd(E_sch)) %>% 
  spread(key = "error_level", value = "sd_E_sch")%>% 
  mutate(E_type = "mean_E_sch")
mean_sd_Eq = results_compiled %>% 
  group_by(Tmax, error_level) %>% 
  summarize(sd_E_q = sd(E_q)) %>% 
  spread(key = "error_level", value = "sd_E_q")%>% 
  mutate(E_type = "mean_E_q")

fit.summary = merge(mean_E, bind_rows(mean_sd_arr, mean_sd_Eq, mean_sd_sch))
colnames(fit.summary)[4:6] = c("error1","error2","error3")

# Basic plotting theme
my_theme = theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Rename variables
fit.summary$E_type[fit.summary$E_type == "mean_E_sch"] = "Schoolfield"
fit.summary$E_type[fit.summary$E_type == "mean_E_arr"] = "Arrhenius"
fit.summary$E_type[fit.summary$E_type == "mean_E_q"] = "Pawar"

# Separate into suboptimal and supraoptimal regions
fit.summary$region = "sub"
fit.summary$region[fit.summary$E_type == "Pawar"] = "supra"
fit.summary$region[fit.summary$Tmax < 32] = "sub"

# create a new point interpolated between 31 and 33 for the Pawar Eq value
# This is just for visualization of Eq
fit.summary %>% 
  subset(E_type == "Pawar") %>% 
  subset(Tmax == 31 | Tmax == 33) %>% 
  group_by(E_type) %>% 
  summarise(Tmax = mean(Tmax),
            E = mean(E),
            error1 = mean(error1),
            error2 = mean(error2),
            error3 = mean(error3)) -> interp 

interp = bind_rows(interp, interp) %>% mutate(region = c("sub", "supra"))
# Make two versions of it, one called sub and one supra
fit.summary = bind_rows(fit.summary, interp)
fit.summary = as.data.frame(fit.summary)

# Make main plot
p1 = ggplot() +
  geom_ribbon(data = fit.summary,aes(x = Tmax, ymin=E-error3, ymax=E+error3, fill = E_type), color = NA, alpha = 0.1) +
  geom_ribbon(data = fit.summary,aes(x = Tmax,ymin=E-error2, ymax=E+error2, fill = E_type), color = NA, alpha = 0.2) +
  geom_ribbon(data = fit.summary,aes(x = Tmax,ymin=E-error1, ymax=E+error1, fill = E_type), color=NA,alpha = 0.3) +
  geom_line(data = fit.summary %>% subset(region == "supra"),aes(x = Tmax, y = E, color = E_type),lty=2) +
  geom_line(data = fit.summary %>% subset(region == "sub"), aes(x = Tmax, y = E, color = E_type),lty = 1) +
  my_theme +
  theme(legend.position = c(0.83,0.83)) +
  theme(legend.title = element_blank()) +
  theme(
    legend.background = element_rect(fill='transparent'), #transparent legend bg
  ) +
  geom_hline(yintercept = 0.66, lty=2) +
  geom_vline(xintercept = 32) +
  ylab(expression(italic(E[a])~"(eV)")) +
  xlab(expression(Temperature~range*","~italic(T[max]-T[min])~"(°C)")) +
  ylim(0,1)

# Make tPC for inset
p2 = ggplotGrob(
  ggplot(tpc, aes(x=temperature,y=performance)) +
  geom_line() +
  xlab("Temperature (°C)") +
  ylab("Performance") +
    my_theme +
    theme(axis.ticks = element_blank(),
          axis.title = element_text(size=8),
          axis.text = element_text(size=6))+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

)

# Place inset
p = p1 + 
  annotation_custom(grob = p2, xmin = 2.5, xmax = 22, ymin = -0.07, ymax = 0.4)

# Write to file
svg("fig5.svg", width = 4, height = 3.25)
grid.arrange(p1)
dev.off()
