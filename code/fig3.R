# Description ----
# Figure 3 plot and analyses from
# Michaletz ST & Garen JC. 2024. Hotter is not (always) better: Embracing unimodal scaling of biological rates with temperature.

# Plot of MST predictions for cross-over of photosynthesis (E = 0.32 eV) and respiration (E = 0.65 eV)

# Initialize ----
#--Libraries
library(ggplot2)
library(ggpubr)
# Uncomment the following lines to install required packages:
# install.packages("ggplot2")
# install.packages("ggpubr")

# create data for the lines
x <- seq(-2, 12, length.out = 100)
y1 <- 0.65 * (x - 5)
y2 <- 0.32 * (x - 5)

# create a data frame with the x and y values for the two lines
df <- data.frame(x, y1, y2)

# Figure 3----
## Panel a: Current MTE prediction----
pan_a <- ggplot(df, aes(x = x)) +
  geom_ribbon(data = df[df$x <= 5, ], aes(x = x, ymin = y1, ymax = y2), fill = "#00B050", alpha = 0.25) +
  geom_ribbon(data = df[df$x >= 5, ], aes(x = x, ymin = y1, ymax = y2), fill = "#ED7D31", alpha = 0.25) +
  geom_line(aes(y = y1), color = "#ED7D31") +
  geom_line(aes(y = y2), color = "#00B050") +
  annotate("text", x = -1.75, y = -2.55, label = "Carbon", color = "#00B050", size = 4, hjust = 0) +
  annotate("text", x = -1.75, y = -3.05, label = "sink", color = "#00B050", size = 4, hjust = 0) +
  annotate("text", x = 10.25, y = 3.1, label = "Carbon", color = "#ED7D31", size = 4, hjust = 0) +
  annotate("text", x = 10.25, y = 2.6, label = "source", color = "#ED7D31", size = 4, hjust = 0) +
  annotate("text", x = 0, y = -0.7, label = expression(italic(E[P]) == 0.32~eV), color = "#00B050", size = 4) +
  annotate("text", x = 8.5, y = 3.85, label = expression(italic(E[R]) == 0.65~eV), color = "#ED7D31", size = 4) +
  xlab(expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste("ln(Photosynthesis or respiration rate)"))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ .*1 , name = expression(paste('Temperature'~(degree*C))))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
pan_a

## Panel b: Revised MTE prediction----
pan_b <- ggplot(df, aes(x = x)) +
  geom_ribbon(data = df, aes(x = x, ymin = y1-1.5, ymax = y1+1.5), fill = "#00B050", alpha = 0.25) +
  geom_line(aes(y = y1-1.5), color = "#ED7D31") +
  geom_line(aes(y = y1+1.5), color = "#00B050") +
  annotate("text", x = 5, y = 0.8, label = "Carbon", color = "#00B050", size = 4, hjust = 0) +
  annotate("text", x = 5, y = 0.3, label = "sink", color = "#00B050", size = 4, hjust = 0) +
  annotate("text", x = 0, y = -0.1, label = expression(italic(E[P]) == 0.65~eV), color = "#00B050", size = 4) +
  annotate("text", x = 9.75, y = -0.2, label = expression(italic(E[R]) == 0.65~eV), color = "#ED7D31", size = 4) +
  xlab(expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste("ln(Photosynthesis or respiration rate)"))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ .*1 , name = expression(paste('Temperature'~(degree*C))))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
pan_b
ggarrange(pan_a, pan_b, common.legend = F, legend = "right", labels = "auto")

## Plot panels together ----
png("Figure_3.png", width = 8,height = 4, res = 300, units = "in")
ggarrange(pan_a, pan_b, common.legend = F, legend = "right", labels = "auto")
dev.off()


