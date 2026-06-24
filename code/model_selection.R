# `survminer` for Kaplan-Meier visualisation
# `survHE` for parametric modelling of survival data
# `muhaz` for non-parametric smoothed hazard estimation
# `landest` for survival estimation from Kaplan-Meier
# `RColorBrewer` for colour mapping
# `pracma` for numerical estimation of restricted mean survival time
# these also load dependencies: ggplot2, ggpubr, survival, flexsurv, dplyr
library(pacman)
p_load(survminer, survHE, muhaz, RColorBrewer, landest, pracma) 

# Load utility functions
source('code/utils.R')

# Save option for AIC-BIC tables
SAVE <- T

# 0. Load digitised data
load('data/digitised_data.Rdata')

# 0.1. Plot KM of internal data (KEYNOTE-006 IA2)
OS.int <- rbind(OS.Pem, OS.Ipi)
OS.int$Treatment <- factor(OS.int$Treatment, levels = c('Pembrolizumab', 'Ipilimumab'))

km_int <- survfit(Surv(Time, Event) ~ Treatment, data = OS.int)

p_km_int <- ggsurvplot(
  km_int,
  data = OS.int,
  censor.shape = '',
  break.time.by = 2,
  risk.table = TRUE,
  xlab = 'Time (months)',
  ylab = 'Overall Survival',
  title = 'Internal data: KEYNOTE-006 IA2 (2015)',
  legend.title = 'Treatment',
  legend.labs = c('Pembrolizumab', 'Ipilimumab')
)

print(p_km_int)

# 0.2. Plot KM of external data (Schadendorf ipilimumab-treatment naive population)
km_ext <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Scha)

p_km_ext <- ggsurvplot(
  km_ext, 
  data = OS.Scha, 
  censor.shape = '', 
  break.time.by = 12, 
  risk.table = TRUE,
  xlab = 'Time (months)', 
  ylab = 'Overall Survival', 
  title = 'External data: Schadendorf treatment naive population (2015)',
  legend = 'none', 
  conf.int = F)

print(p_km_ext)

# 1. Internal model for pembrolizumab arm
# Non-parametric smoothed hazard - increasing then decreasing shape
haz_Pem <- muhaz(OS.Pem$Time, OS.Pem$Event)

plot(haz_Pem, xlab = 'Time (months)', main = 'Pembrolizumab - Smoothed Hazard')

# 1.1. Standard parametric model
formula <- Surv(Time, Event) ~ 1
mods <- c('exp', 'weibull', 'gompertz', 'gengamma', 'loglogistic', 'lognormal')

m_Pem_param <- fit.models(formula = formula, data = OS.Pem, distr = mods)

models <- list(
  'Generalised Gamma' = m_Pem_param$models$`Gen. Gamma`,
  'Log-Logistic' = m_Pem_param$models$`log-Logistic`,
  'Log-Normal' = m_Pem_param$models$`log-Normal`
)

# Visual inspection of hazard plot
plot_hazard(
  models_to_plot = models, 
  haz = haz_Pem, 
  title = 'Pembrolizumab - Standard Parametric Models')

# Visual inspection of survival plot
km_Pem <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Pem)

plot_survival(
  models_to_plot = models, 
  km = km_Pem, 
  time_points = seq(0,20,0.1), 
  title = 'Pembrolizumab - Standard Parametric Models')

# AIC, BIC
aic_bic_summary_Pem_param <- data.frame(
  AIC = sapply(m_Pem_param$models, AIC), 
  BIC = sapply(m_Pem_param$models, BIC))

print(aic_bic_summary_Pem_param)

# 1.2. 1-knot cubic spline model
m_Pem_spline_1_hazard <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'hazard', 
  k = 1)

m_Pem_spline_1_odds <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'odds', 
  k = 1)

# Error in optim using 'BFGS', thus changing to 'Nelder-Mead'
m_Pem_spline_1_normal <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'normal', 
  k = 1, 
  method = "Nelder-Mead")

models <- list(
  '1-knot spline hazard' = m_Pem_spline_1_hazard,
  '1-knot spline odds' = m_Pem_spline_1_odds,
  '1-knot spline normal' = m_Pem_spline_1_normal
)

# Visual inspection of hazard
plot_hazard(
  models_to_plot = models, 
  haz = haz_Pem, 
  title = 'Pembrolizumab - 1-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(
  models_to_plot = models, 
  km = km_Pem, 
  time_points = seq(0,20,0.1), 
  title = 'Pembrolizumab - 1-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Pem_spline_1 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_1)

# 1.3. 2-knot cubic spline model
m_Pem_spline_2_hazard <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'hazard', 
  k = 2)

m_Pem_spline_2_odds <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'odds', 
  k = 2)

m_Pem_spline_2_normal <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'normal', 
  k = 2)

models <- list(
  '2-knot spline hazard' = m_Pem_spline_2_hazard,
  '2-knot spline odds' = m_Pem_spline_2_odds,
  '2-knot spline normal' = m_Pem_spline_2_normal
)

# Visual inspection of hazard
plot_hazard(
  models_to_plot = models, 
  haz = haz_Pem, 
  title = 'Pembrolizumab - 2-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(
  models_to_plot = models, 
  km = km_Pem, 
  time_points = seq(0,20,0.1), 
  title = 'Pembrolizumab - 2-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Pem_spline_2 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_2)

# 1.4. 3-knot cubic spline model
m_Pem_spline_3_hazard <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'hazard', 
  k = 3)

m_Pem_spline_3_odds <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'odds', 
  k = 3)

m_Pem_spline_3_normal <- flexsurvspline(
  formula = formula, 
  data = OS.Pem, 
  scale = 'normal', 
  k = 3)

models <- list(
  '3-knot spline hazard' = m_Pem_spline_3_hazard,
  '3-knot spline odds' = m_Pem_spline_3_odds,
  '3-knot spline normal' = m_Pem_spline_3_normal
)

# Visual inspection of hazard
plot_hazard(
  models_to_plot = models, 
  haz = haz_Pem, 
  title = 'Pembrolizumab - 3-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(
  models_to_plot = models, 
  km = km_Pem, 
  time_points = seq(0, 20, 0.1), 
  title = 'Pembrolizumab - 3-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Pem_spline_3 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_3)

# 2. Internal model for Ipilimumab arm 
# Non-parametric smoothed hazard - increasing then decreasing
haz_Ipi <- muhaz(OS.Ipi$Time, OS.Ipi$Event, bw.smooth = 3)

plot(haz_Ipi, xlab='Time (months)', main='Ipilimumab - Smoothed Hazard')

# 2.1. Standard parametric model
m_Ipi_param <- fit.models(formula = formula, data = OS.Ipi, distr = mods)

models <- list(
  'Generalised Gamma' = m_Ipi_param$models$`Gen. Gamma`,
  'Log-Logistic' = m_Ipi_param$models$`log-Logistic`,
  'Log-Normal' = m_Ipi_param$models$`log-Normal`
)

# Visual inspection of hazard
plot_hazard(
  models_to_plot = models, 
  haz = haz_Ipi, 
  title = 'Ipilimumab - Standard Parametric Models')

# Visual inspection of survival
km_Ipi <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Ipi)

plot_survival(
  models_to_plot = models, 
  km = km_Ipi, 
  time_points = seq(0,20,0.1), 
  title = 'Ipilimumab - Standard Parametric Models')

# AIC, BIC
aic_bic_summary_Ipi_param <- data.frame(
  AIC = sapply(m_Ipi_param$models, AIC), 
  BIC = sapply(m_Ipi_param$models, BIC))

print(aic_bic_summary_Ipi_param)

# 2.2 1-knot spline model
m_Ipi_spline_1_hazard <- flexsurvspline(
  formula = formula, 
  data = OS.Ipi, 
  scale = 'hazard', 
  k = 1)

m_Ipi_spline_1_odds <- flexsurvspline(
  formula = formula, 
  data = OS.Ipi, 
  scale = 'odds', 
  k = 1)

m_Ipi_spline_1_normal <- flexsurvspline(
  formula = formula, 
  data = OS.Ipi, 
  scale = 'normal', 
  k = 1)

models <- list(
  '1-knot spline hazard' = m_Ipi_spline_1_hazard,
  '1-knot spline odds' = m_Ipi_spline_1_odds,
  '1-knot spline normal' = m_Ipi_spline_1_normal
)

# Visual inspection of hazard
plot_hazard(
  models_to_plot = models, 
  haz = haz_Ipi, 
  title = 'Ipilimumab - 1-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(
  models_to_plot = models, 
  km = km_Ipi, 
  time_points = seq(0, 20, 0.1), 
  title = 'Ipilimumab - 1-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Ipi_spline_1 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Ipi_spline_1)

# 3. External model for both arms
haz_Scha <- muhaz(OS.Scha$Time, OS.Scha$Event, bw.smooth = 6)

plot(haz_Scha, xlab='Time (months)', main='Schadendorf Study (2015) -  Smoothed Hazard')
abline(v = 13.84, col = 'grey', lty = 2)

# 3.1. Rebased standard parametric model
rebased_time <- 13.84 # rebased_time = median follow-up of KEYNOTE006 IA2

OS.Scha.rebased <- subset(OS.Scha, Time > rebased_time)
OS.Scha.rebased$Time <-  OS.Scha.rebased$Time - rebased_time

m_Scha_rebased = fit.models(formula = formula, data = OS.Scha.rebased, distr = mods)

models <- list(
  'Rebased Weibull' = m_Scha_rebased$models$`Weibull (AFT)`,
  'Rebased Gompertz' = m_Scha_rebased$models$Gompertz,
  'Rebased Gen Gamma' = m_Scha_rebased$models$`Gen. Gamma`,
  'Rebased Loglogistic' = m_Scha_rebased$models$`log-Logistic`,
  'Rebased Lognormal' = m_Scha_rebased$models$`log-Normal`
)

# Visual inspection of hazard
plot_rebased_hazard(
  models_to_plot = models, 
  haz = haz_Scha, 
  rebased_time = rebased_time, 
  title = 'Schadendorf (2015) - Hazard of Parametric Models Rebased at 13.84 Months')

# Visual inspection of survival
km_Scha <- survfit(Surv(Time, Event) ~ 1, data = OS.Scha)

plot_rebased_survival(
  models_to_plot = models, 
  km = km_Scha, 
  rebased_time = rebased_time, 
  time_points = seq(0, 84, by = 0.1), 
  title = 'Schadendorf (2015) - Survival of Parametric Models Rebased at 13.84 Months')

# AIC, BIC
aic_bic_summary_Scha_rebased <- data.frame(
  AIC = sapply(m_Scha_rebased$models, AIC), 
  BIC = sapply(m_Scha_rebased$models, BIC))

print(aic_bic_summary_Scha_rebased)

# Save AIC-BIC tables
if (SAVE == T) {
  write.csv(
    rbind(
      aic_bic_summary_Pem_param, 
      aic_bic_summary_Pem_spline_1, 
      aic_bic_summary_Pem_spline_2, 
      aic_bic_summary_Pem_spline_3
    ), 
    'tables/model selection/aic_bic_int_pem_selection.csv', 
    row.names = T)
  write.csv(
    rbind(
      aic_bic_summary_Ipi_param, 
      aic_bic_summary_Ipi_spline_1
    ), 
    'tables/model selection/aic_bic_int_ipi_selection.csv', 
    row.names = T)
  write.csv(
    aic_bic_summary_Scha_rebased, 
    'tables/model selection/aic_bic_ext_selection.csv', 
    row.names = T)
}