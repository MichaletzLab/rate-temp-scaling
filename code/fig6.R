# Code to generate Figure 6 - Arroyo et al. model parameter estimates
# Sensitivity to noise and T range
# Michaletz ST & Garen JC. 2024. Hotter is not (always) better: Embracing unimodal scaling of biological rates with temperature.

library(rTPC)
library(tidyverse)
library(ggplot2)
library(nls.multstart)
library(cowplot)

# Basic plotting theme
my_theme = theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Model from Arroyo et al.
arroyo_model = function(TC, dH, dC, lnB0) {
  
  # TC = temperature (C)
  # dH = change in enthalphy (J/mol)
  # dC = change in heat capacity (J/K mol)
  # lnB0 = normalization constant

  R = 8.314462618 # universal gas constant (J/K mol)
  
  # Convert to K
  TK = TC + 273.15

  ln.rate = lnB0 - (dH/R)*(1/TK) - (dC/R)*log(1/TK)
  return(exp(ln.rate))
  
}


# Loop over T domain width
results = data.frame()

tmin_abs = 5
tmax_abs = 50
incr = 0.5
trange_abs = seq(tmin_abs, tmax_abs, incr)

for(i in trange_abs) {
  
  # Loop over T domain location
  for(j in seq(0,50-i,incr)) {
    # Generate TPC
    tmin = j
    tmax = j+i
    trange = seq(tmin, tmax, incr)
    cur.tpc = data.frame(temperature = trange,
                         performance = pawar_2018(trange, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10),
                         T_ref = 25)
 
    # Fit Arroyo model
    fit = nls_multstart(performance ~ arroyo_model(TC = temperature, dH, dC, lnB0),
                        data = cur.tpc,
                        iter = 1000,
                        start_lower = c(dH = 10000, dC = -1000, lnB0 = -100),
                        start_upper = c(dH = 100000, dC = 1000, lnB0 = 100),
                        supp_errors = 'Y',
                        na.action = na.omit)
    print(coef(fit)[3])
    
    # Compile results
    cur_results = data.frame(Tmin = j,
                             Tmax = j+i,
                             Width = i,
                             dH = coef(fit)[1],
                             dC = coef(fit)[2])
    results = bind_rows(results, cur_results)
    
  }
}

results$Mean_T = (results$Tmin+results$Tmax)/2
results = subset(results, Mean_T%%(incr) == 0)


# Plot heatmap for dH
p1 = ggplot(results, aes(x=Width,y=Mean_T,fill=dH/1000, color=dH/1000)) + 
  geom_tile() + 
  my_theme +
  xlab(expression(Temperature~range*","~italic(T[max]-T[min])~"(°C)")) +
  ylab(expression("Temperature range-location, <"~italic(T)~"> (°C)")) +
  theme(legend.position = c(0.67,0.86),
        legend.direction = "horizontal") +
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),  
        legend.text=element_text(size=rel(0.7), angle=45, hjust=1),
        legend.title = element_text(size=rel(0.9), vjust = 0.9),
        legend.key.size = unit(0.5,"cm")) +
  #theme(legend.title = element_text("hello")) +
  scale_colour_distiller(name="ΔH (kJ)",palette = "Spectral", direction = 1)+
  scale_fill_distiller(name="ΔH (kJ)",palette = "Spectral", direction = 1) +
  labs(tag = "A")


# Plot heatmap for dC
p2 = ggplot(results, aes(x=Width,y=Mean_T,fill=dC/1000, color=dC/1000)) + 
  geom_tile() + 
  my_theme +
  xlab(expression(Temperature~range*","~italic(T[max]-T[min])~"(°C)")) +
  ylab(expression("Temperature range-location, <"~italic(T)~"> (°C)")) +
  theme(legend.position = c(0.65,0.89),
        legend.direction = "horizontal") +
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=rel(0.9), vjust = 0.9),
        legend.key.size = unit(0.5,"cm")) +
  scale_colour_distiller(name=expression("ΔC (kJ K"^-1*")"),palette = "Spectral", direction = 1)+
  scale_fill_distiller(name=expression("ΔC (kJ K"^-1*")"),palette = "Spectral", direction = 1) +
  labs(tag="B")

## Part 2

# Function for parallel processing curve fitting
fit_curves_parallel = function(dat) {
  
  # Fit curve with NLS multstart
  fit = nls_multstart(Photo ~ arroyo_model(TC = Tleaf, dH, dC, lnB0),
                      data = dat,
                      iter = 1000,
                      start_lower = c(dH = 10000, dC = -1000, lnB0 = -100),
                      start_upper = c(dH = 100000, dC = 1000, lnB0 = 100),
                      supp_errors = 'N',
                      na.action = na.omit)
  
  # If the fit was successful, extract parameters estimates and SE
  if(typeof(fit) != "NULL") { 
    dH <- coef(fit)["dH"]
    dH_SE <- summary(fit)$coefficients[,'Std. Error']['dH']
    dC <- coef(fit)["dC"]
    dC_SE <- summary(fit)$coefficients[,'Std. Error']['dC']
    AIC <- AIC(fit)
    r_sq = 1-sum(resid(fit)^2)/sum((dat$Photo-mean(dat$Photo))^2)
    
    # Otherwise, likely a convergence failure. All parameters set to NA
  } else { 
    dH = dH_SE = dC = dC_SE = AIC = r_sq = NA
    dat$failure_status = "Convergence failure"
  }
  
  # Make data frame with results and return
  results = data.frame(select(dat[1,],-Tleaf, -Photo), dH, dH_SE, 
                       dC, dC_SE, AIC, r_sq)
  return(results)
  
}

# Generate TPC to use
tr = seq(from=10,to=40,length.out=10)
tpc.curve = data.frame(temperature = tr,
                       performance = pawar_2018(tr, r_tref = 1, e = 0.66, eh = 1.15, topt = 25.3, tref = 10),
                       T_ref = 25)

## Introduce gaussian noise at progressively larger levels and fit curves
nrep = 100 # Number of repetitions at each noise level
noise.levels = seq(0.001,0.03, 0.001) # SD of gaussian noise
ntot = nrep*length(noise.levels)

len = dim(tpc.curve)[1]

# copy each curve a number of times equal to noise levels times nrep
tpc.curves = tpc.curve[rep(1:len,ntot),]
tpc.curves$noise.level = rep(noise.levels, each=(nrep*len))
tpc.curves$curveID = rep(1:ntot, each=len)

# Function for adding noise to TPCs
func1 <- function(x, y) {rnorm(1, mean = x, sd = y) }
list1 <- Map(func1, 0, tpc.curves$noise.level)

# Apply noise
tpc.curves$A_noise = tpc.curves$performance + as.numeric(Map(func1, 0, tpc.curves$noise.level))

# Bind tpc without noise added
tpc.curves = bind_rows(tpc.curve %>% mutate(A_noise = performance, noise.level = 0, curveID = 1), 
                       tpc.curve %>% mutate(A_noise = performance, noise.level = 0, curveID = 2),
                       tpc.curves)

tpc.dat = tpc.curves %>% select(curveID, Tleaf = temperature, Photo = A_noise, noise.level) %>% mutate(type = "High")

# Compile tpcs and apply curveIDs
dat = tpc.dat 
dat$curveID = dat %>% group_by(type, curveID) %>% group_indices()

## Fit curves with parallel processing function
library(parallel)
print("Fitting curves, please wait...")
dat2 = list()
j = 1
for (i in unique(dat$curveID)) {
  cur = subset(dat, curveID == i)
  dat2[[j]] = as.data.frame(cur)
  j = j + 1
  
}
system.time({
  numcores = detectCores()
  clust <- makeCluster(numcores)
  clusterExport(clust, "nls_multstart")
  clusterExport(clust, "arroyo_model")
  clusterExport(clust, "select")
  results2 = parLapply(clust, dat2, fit_curves_parallel)
})
results2 = as.data.frame(do.call(rbind, results2))

print("Curve fitting complete")

# Collect results and compute SD
results.summary = results2 %>% 
  group_by(noise.level, type) %>% 
  summarize(mean_dH = mean(dH), sd_dH = sd(dH),
            mean_dC = mean(dC), sd_dC = sd(dC))

# Reshape data
results.plot.data = results.summary %>% 
  gather(key = "parameter", value = "value", mean_dH, 
         mean_dC)

# Rearrange dataframe
results.sd.data = results.plot.data
results.sd.data$sd = NA
results.sd.data[results.sd.data$parameter == "mean_dH",]$sd = results.sd.data[results.sd.data$parameter == "mean_dH",]$sd_dH
results.sd.data[results.sd.data$parameter == "mean_dC",]$sd = results.sd.data[results.sd.data$parameter == "mean_dC",]$sd_dC

# Make shaded interval regions
results.sd.data = results.sd.data %>% mutate(upper = value+sd, lower = value-sd) 

# Fix labels
plot.data.2 = results.plot.data %>% rename(Density = type)
sd.data.2 = results.sd.data %>% rename(Density = type)

# Panel 3: dH noise response
a = subset(plot.data.2, parameter == "mean_dH")
b = subset(sd.data.2, parameter == "mean_dH")
p3 = ggplot(a, aes(x=100*noise.level, y = value/1000, color = Density)) +
  geom_line() +  
  geom_ribbon(data=b, aes(x =100*noise.level, ymin = lower/1000, ymax = upper/1000, fill=Density), color=NA,alpha = 0.2)+
  my_theme +
  xlab("Noise level (%)") +
  ylab("ΔH (kJ)") +
  labs(tag = "C")

# Panel 4: dC noise respones
a = subset(plot.data.2, parameter == "mean_dC")
b = subset(sd.data.2, parameter == "mean_dC")
p4 = ggplot(a, aes(x=100*noise.level, y = value/1000, color = Density)) +
  geom_line() +  
  geom_ribbon(data=b, aes(x =100*noise.level, ymin = lower/1000, ymax = upper/1000, fill=Density), color=NA,alpha = 0.2)+
  my_theme +
  xlab("Noise level (%)") +
  ylab(expression("ΔC (kJ K"^-1*")")) +
  labs(tag = "D")

# Save plot to files
svg("Figure_6.svg",width=7,height=6)
plot_grid(p1,p2,p3,p4,
          align="hv")
dev.off()
postscript("Figure_6.eps", width = 7, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
plot_grid(p1,p2,p3,p4,
          align="hv")
dev.off()

