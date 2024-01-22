# Figure S1 plot and analyses from
# Michaletz ST & Garen JC. 2024. Hotter is not (always) better: Embracing unimodal scaling of biological rates with temperature.

library(rTPC)
library(tidyverse)
library(ggplot2)

# Basic plotting theme
my_theme = theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Generate a full TPC
tpc = data.frame(temperature = 0:50,
                 performance = pawar_2018(0:50, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10))


# Initiate empty holder variable
fitting_data = c()

# Loop over temperature ranges from 5-45
for (i in seq(5,45,2)) {
  
  # Generate 10 points equally spaced over the temperature range
  cur = data.frame(temperature = i*(0:9)/9,
                   performance = pawar_2018( i*(0:9)/9, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10),
                   t_range = i)
  
  # Bind together
  fitting_data = bind_rows(fitting_data, cur)
  
}


# Build plot
png("fig_S1.png", height = 6, width = 6, units="in", res=300)
ggplot(fitting_data, aes(x = temperature, y = performance)) +
  geom_line(data = tpc %>% subset(temperature <= 45), color = "grey", linewidth=1.5) +
  geom_point() + 
  my_theme +
  #theme(strip.placement) +
  xlab("Temperature (°C)") +
  ylab("Biological rate") +
  facet_wrap(~t_range,ncol=4, strip.position="right",
             labeller = as_labeller(c("5" = "T[range]==5*~degree*C",
                                      "7" = "T[range]==7*~degree*C",
                                      "9" = "T[range]==9*~degree*C",
                                      "11" = "T[range]==11*~degree*C",
                                      "13" = "T[range]==13*~degree*C",
                                      "15" = "T[range]==15*~degree*C",
                                      "17" = "T[range]==17*~degree*C",
                                      "19" = "T[range]==19*~degree*C",
                                      "21" = "T[range]==21*~degree*C",
                                      "23" = "T[range]==23*~degree*C",
                                      "25" = "T[range]==25*~degree*C",
                                      "27" = "T[range]==27*~degree*C",
                                      "29" = "T[range]==29*~degree*C",
                                      "31" = "T[range]==31*~degree*C",
                                      "33" = "T[range]==33*~degree*C",
                                      "35" = "T[range]==35*~degree*C",
                                      "37" = "T[range]==37*~degree*C",
                                      "39" = "T[range]==39*~degree*C",
                                      "41" = "T[range]==41*~degree*C",
                                      "43" = "T[range]==43*~degree*C",
                                      "45" = "T[range]==45*~degree*C"),
                                    default = label_parsed))
dev.off()
