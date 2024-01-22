# Code for Fig. 4 - heatmap of Arrhenius Ea values

library(rTPC)
library(tidyverse)
library(ggplot2)
library(nls.multstart)
library(gridExtra)


# Basic plotting theme
my_theme = theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


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

# Initialize results dataframe
results = data.frame()

# Set minimum and maximum T values and increment size
tmin_abs = 5
tmax_abs = 50
incr = 0.5
trange_abs = seq(tmin_abs, tmax_abs, incr)

# Loop over T domain width
for(i in trange_abs) {
  
  # Loop over T domain location
  for(j in seq(0,50-i,incr)) {
    
    # Generate TPC w/ SS model over current T range and range-location
    tmin = j
    tmax = j+i
    trange = seq(tmin, tmax, incr)
    cur.tpc = data.frame(temperature = trange,
                         performance = pawar_2018(trange, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10))
    
    # Fit Arrhenius model to TPC
    fit = nls(performance ~ arrhenius(temperature, J_ref, E, T_ref = 10),
              data = cur.tpc,
              start = list(J_ref = 1, E = 0.4))
    
    # Save results
    cur_results = data.frame(Tmin = j,
                             Tmax = j+i,
                             Width = i,
                             J_ref = coef(fit)[1],
                             E_arr = coef(fit)[2])
    results = bind_rows(results, cur_results)
  }
}

# Compute mean T
results$Mean_T = (results$Tmin+results$Tmax)/2
results = subset(results, Mean_T%%(incr) == 0)


# Plot heatmap panel A
p1 = ggplot(results, aes(x=Width,y=Mean_T,fill=E_arr, color=E_arr)) + 
  geom_tile() + 
  my_theme +
  xlab(expression(Temperature~range*","~italic(T[max]-T[min])~"(°C)")) +
  ylab(expression("Temperature range-location, <"~italic(T)~"> (°C)")) +
  theme(legend.position = c(0.75,0.85),
        legend.direction = "horizontal") +
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_colour_distiller(name=expression(italic(E[a*","*A])),palette = "Spectral", direction = 1)+
  scale_fill_distiller(name=expression(italic(E[a*","*A])),palette = "Spectral", direction = 1)+
  labs(tag="A")

# Generate TPC
tpc = data.frame(temperature = 0:50,
                 performance = pawar_2018(0:50, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10))

# Panel b showing shape of TPC
p2 = ggplot(tpc, aes(x=temperature,y=performance)) +
    geom_line() +
    xlab("Temperature (°C)") +
    ylab("Biological rate") +
    my_theme +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ) +
  labs(tag = "B")

# Panel C showing density of Ea values
p3 = ggplot(data = results, aes(x = E_arr)) +
  geom_density() +
  my_theme +
  xlab(expression(italic(E[a*","*A])~"(eV)")) +
  ylab("Density") +
  labs(tag = "C")

svg("fig4.svg", width = 7, height = 4)
grid.arrange(p1,p2,p3,ncol=2,
             layout_matrix = cbind(c(1,1), c(1,1), c(2,3)))
dev.off()
