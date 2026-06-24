library(pacman)
p_load(survHE)
source('code/utils.R')
load('data/digitised_data.Rdata')

formula <- Surv(Time, Event) ~ 1
m_Pem <- fit.models(formula = formula, data = OS.Pem.7y, distr = 'gengamma')
m_Ipi <- fit.models(formula = formula, data = OS.Ipi.7y, distr = 'gengamma')
m_Scha <- fit.models(formula = formula, data = OS.Scha, distr = 'gompertz')

time_horizon <- 84 
dt = 0.1
t_seq <- seq(0, time_horizon, dt)

h_Pem <- as.data.frame(summary(m_Pem$models$`Gen. Gamma`, type = 'hazard', t = t_seq))
h_Ipi <- as.data.frame(summary(m_Ipi$models$`Gen. Gamma`, type = 'hazard', t = t_seq))
h_Scha <- as.data.frame(summary(m_Scha$models$Gompertz, type = 'hazard', t = t_seq))
h_Pem$est[1] <- 0
h_Ipi$est[1] <- 0

t1 <- 12
t2 <- 60
a <- 7
b <- 3

# hazard
h_Pem_blended <- compute_blended_hazard(h_Pem, h_Scha, t1, t2, a, b, t_seq) 
h_Ipi_blended <- compute_blended_hazard(h_Ipi, h_Scha, t1, t2, a, b, t_seq)

hazplot_all <- ggplot() +
  geom_line(data = h_Pem, aes(x = time, y = est, color = 'Fitted Internal Hazard (Arm 1)'), linewidth = 1, linetype = 'dashed') + 
  geom_line(data = h_Ipi, aes(x = time, y = est, color = 'Fitted Internal Hazard (Arm 0)'), linewidth = 1, linetype = 'dashed') +
  geom_line(data = h_Scha, aes(x = time, y = est, color = 'Fitted External Hazard'), linewidth = 1, linetype = 'dashed') +
  geom_line(data = h_Pem_blended, aes(x = time, y = est, color = 'Blended Hazard (Arm 1)'), linewidth = 1) +
  geom_line(data = h_Ipi_blended, aes(x = time, y = est, color = 'Blended Hazard (Arm 0)'), linewidth = 1) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) + ylim(0, 0.065) +
  scale_color_brewer(palette = 'Set1', name = '') + 
  coord_cartesian(clip = 'off') +
  theme_classic() + theme(axis.ticks.x = element_blank(), 
                          axis.line.x = element_line(), 
                          legend.position = c(0.75, 0.8), 
                          text = element_text(size = 12), 
                          legend.key.width = unit(1, 'cm')) +
  labs(x = 'Time', y = 'Hazard')

hazplot_arm0 <- ggplot() +
  geom_line(data = h_Ipi, aes(x = time, y = est, color = 'Fitted Internal Hazard'), linewidth = 1, linetype = 'dashed') +
  geom_line(data = h_Scha, aes(x = time, y = est, color = 'Fitted External Hazard'), linewidth = 1, linetype = 'dashed') +
  geom_line(data = h_Ipi_blended, aes(x = time, y = est, color = 'Blended Hazard'), linewidth = 1) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) + ylim(0, 0.065) +
  scale_color_manual(
    name = "",
    values = c(
      "Fitted Internal Hazard" = "#E41A1C",
      "Fitted External Hazard" = "#4DAF4A",
      "Blended Hazard" = "black"
    )
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic() + theme(axis.ticks.x = element_blank(), 
                          axis.line.x = element_line(), 
                          legend.position = c(0.75, 0.8), 
                          text = element_text(size = 12), 
                          legend.key.width = unit(1, 'cm')) +
  labs(x = 'Time', y = '')

hazplot_arm1 <- ggplot() +
  geom_line(data = h_Pem, aes(x = time, y = est, color = 'Fitted Internal Hazard'), linewidth = 1, linetype = 'dashed') + 
  geom_line(data = h_Scha, aes(x = time, y = est, color = 'Fitted External Hazard'), linewidth = 1, linetype = 'dashed') +
  geom_line(data = h_Pem_blended, aes(x = time, y = est, color = 'Blended Hazard'), linewidth = 1) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) + ylim(0, 0.065) +
  scale_color_manual(
    name = "",
    values = c(
      "Fitted Internal Hazard" = "#E41A1C",
      "Fitted External Hazard" = "#4DAF4A",
      "Blended Hazard" = "black"
    )
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic() + theme(axis.ticks.x = element_blank(), 
                          axis.line.x = element_line(), 
                          legend.position = c(0.75, 0.8), 
                          text = element_text(size = 12), 
                          legend.key.width = unit(1, 'cm')) +
  labs(x = 'Time', y = '')

# HR
HR_df <- merge(h_Pem_blended, h_Ipi_blended, by = "time", suffixes = c("_arm1", "_arm0"))
HR_df$HR <- HR_df$est_arm1 / HR_df$est_arm0
HR_df_sub <- subset(HR_df, time > 0.2)

HRplot <- ggplot(HR_df_sub, aes(x = time, y = HR)) + 
  geom_line(color = 1, linewidth = 1.2) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0.5, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0.5, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) + ylim(0.5, 1) +
  coord_cartesian(clip = 'off') +
  theme_classic() + theme(axis.ticks.x = element_blank(), 
                          axis.line.x = element_line(), 
                          text = element_text(size = 12)) +
  labs(x = 'Time', y = 'HR')

# survival
compute_survival <- function(df) {
  df <- df[order(df$time), ]
  delta_t <- c(0, diff(df$time))  # time step
  df$cumhaz <- cumsum(df$est * delta_t)
  df$surv <- exp(-df$cumhaz)
  return(df)
}

S_Pem <- compute_survival(h_Pem_blended)
S_Ipi <- compute_survival(h_Ipi_blended)

survplot_all <- ggplot() +
  geom_line(data = S_Pem, aes(x = time, y = surv, color = "Arm 1"), linewidth = 1) +
  geom_line(data = S_Ipi, aes(x = time, y = surv, color = "Arm 0"), linewidth = 1) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) +
  ylim(0, 1) +
  labs(x = "Time", y = "Survival", title = "") +
  scale_color_brewer(palette = 'Set1', name = '') + coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.line.x = element_line(), 
    text = element_text(size = 12),
    legend.position = c(0.8, 0.85),
    legend.key.width = unit(1, 'cm')
  )

survplot_arm0 <- ggplot() +
  geom_line(data = S_Ipi, aes(x = time, y = surv), color = "#E41A1C", linewidth = 1) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) +
  ylim(0, 1) +
  labs(x = "Time", y = "", title = "") +
  scale_color_brewer(palette = 'Set1', name = '') + coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.line.x = element_line(), 
    text = element_text(size = 12),
    legend.position = c(0.8, 0.85),
    legend.key.width = unit(1, 'cm')
  )

survplot_arm1 <- ggplot() +
  geom_line(data = S_Pem, aes(x = time, y = surv), color = "#E41A1C", linewidth = 1) +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  annotate("text", x = t1, y = 0, label = expression(t[1]), vjust = 3, size = 4.5) +
  annotate("text", x = t2, y = 0, label = expression(t[2]), vjust = 3, size = 4.5) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12), labels = rep("", 8)) +
  ylim(0, 1) +
  labs(x = "Time", y = "", title = "") +
  scale_color_brewer(palette = 'Set1', name = '') + coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.line.x = element_line(), 
    text = element_text(size = 12),
    legend.position = c(0.8, 0.85),
    legend.key.width = unit(1, 'cm')
  )

library(patchwork)
combined_plot_1 <- hazplot_arm0 + hazplot_arm1 + HRplot + plot_layout(ncol = 3)
combined_plot_2 <- hazplot_arm0 + hazplot_arm1 + HRplot + survplot_all + plot_layout(ncol = 2, nrow = 2)
combined_plot_2

############
library(gridExtra)
library(grid)

# Column names
col1 <- textGrob("Arm 0", gp=gpar(fontsize=12, fontface="bold"))
col2 <- textGrob("Arm 1", gp=gpar(fontsize=12, fontface="bold"))

# Row names
row1 <- textGrob("Hazard", rot = 90, gp=gpar(fontsize=12, fontface="bold"))
row2 <- textGrob("Survival", rot = 90, gp=gpar(fontsize=12, fontface="bold"))

grid.arrange(
  arrangeGrob(nullGrob(), col1, col2, ncol = 3, widths = c(0.5, 6, 6)),
  arrangeGrob(row1, hazplot_arm0, hazplot_arm1, ncol = 3, widths = c(0.5, 6, 6)),
  arrangeGrob(row2, survplot_arm0, survplot_arm1, ncol = 3, widths = c(0.5, 6, 6)),
  nrow = 3,
  heights = c(0.5, 6, 6)
)

HRplot

