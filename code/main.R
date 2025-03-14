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

# 0. Load digitised data
load('data/digitised_data.Rdata')

# 0.1. Plot KM of internal data (KEYNOTE-006 IA2)
OS.int <- rbind(OS.Pem, OS.Ipi)
OS.int$Treatment <- factor(OS.int$Treatment, levels = c('Pembrolizumab', 'Ipilimumab'))
km_int <- survfit(Surv(Time, Event) ~ Treatment, data = OS.int)

p_km_int <- ggsurvplot(km_int, data = OS.int, 
                       censor.shape = '', break.time.by = 2, risk.table = TRUE,
                       xlab = 'Time (months)', ylab = 'Overall Survival', 
                       legend.title = 'Treatment', legend.labs = c('Pembrolizumab', 'Ipilimumab'))
p_km_int$plot <-  p_km_int$plot + 
  labs(title = 'Internal data: KEYNOTE-006 IA2 (2015)')

print(p_km_int)

# 0.2. Plot KM of external data (Schadendorf ipilimumab-treatment naive population)
km_ext <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Scha)

p_km_ext <- ggsurvplot(km_ext, data = OS.Scha, 
                       censor.shape = '', break.time.by = 12, risk.table = TRUE,
                       xlab = 'Time (months)', ylab = 'Overall Survival', 
                       legend = 'none', conf.int = F)
p_km_ext$plot <-  p_km_ext$plot + 
  labs(title = 'External data: Schadendorf treatment naive population (2015)')

print(p_km_ext)

# 1. Internal model for pembrolizumab arm
# Non-parametric smoothed hazard - increasing then decreasing shape
haz_Pem <- muhaz(OS.Pem$Time, OS.Pem$Event)
plot(haz_Pem, xlab = 'Time (months)', main = 'Pembrolizumab - Smoothed Hazard')

# 1.1. Standard parametric model
formula <- Surv(Time,Event) ~ 1
mods <- c('exp', 'weibull', 'gompertz', 'gengamma', 'loglogistic', 'lognormal')
m_Pem_param <- fit.models(formula = formula, data = OS.Pem, distr = mods)

# AIC, BIC
aic_bic_summary_Pem_param <- data.frame(AIC = sapply(m_Pem_param$models, AIC), BIC = sapply(m_Pem_param$models, BIC))
print(aic_bic_summary_Pem_param)

# Visual inspection of hazard plot

models <- list(
  'Generalised Gamma' = m_Pem_param$models$`Gen. Gamma`,
  'Log-Logistic' = m_Pem_param$models$`log-Logistic`,
  'Log-Normal' = m_Pem_param$models$`log-Normal`
)

plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - Standard Parametric Models')

# Visual inspection of survival plot
km_Pem <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Pem)
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,24,0.1), title = 'Pembrolizumab - Standard Parametric Models')

# 1.2. 1-knot cubic spline model
m_Pem_spline_1_hazard <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'hazard', k = 1)
m_Pem_spline_1_odds <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'odds', k = 1)

# AIC, BIC
models <- list(
  '1-knot spline hazard' = m_Pem_spline_1_hazard,
  '1-knot spline odds' = m_Pem_spline_1_odds
)

aic_bic_summary_Pem_spline_1 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(Model = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_1)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - 1-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,24,0.1), title = 'Pembrolizumab - 1-Knot Cubic Spline Models')

# 1.3. 2-knot cubic spline model
m_Pem_spline_2_hazard <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'hazard', k = 2)
m_Pem_spline_2_odds <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'odds', k = 2)
m_Pem_spline_2_normal <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'normal', k = 2)

# AIC, BIC
models <- list(
  '2-knot spline hazard' = m_Pem_spline_2_hazard,
  '2-knot spline odds' = m_Pem_spline_2_odds,
  '2-knot spline normal' = m_Pem_spline_2_normal
)

aic_bic_summary_Pem_spline_2 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(Model = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_2)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - 2-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,24,0.1), title = 'Pembrolizumab - 2-Knot Cubic Spline Models')

# 1.4. 3-knot cubic spline model
m_Pem_spline_3_hazard <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'hazard', k = 3)
m_Pem_spline_3_odds <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'odds', k = 3)
m_Pem_spline_3_normal <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'normal', k = 3)

# AIC, BIC
models <- list(
  '3-knot spline hazard' = m_Pem_spline_3_hazard,
  '3-knot spline odds' = m_Pem_spline_3_odds,
  '3-knot spline normal' = m_Pem_spline_3_normal
)

aic_bic_summary_Pem_spline_3 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(Model = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_3)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - 3-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,24,0.1), title = 'Pembrolizumab - 3-Knot Cubic Spline Models')

# 2. Internal model for Ipilimumab arm
# Non-parametric smoothed hazard - increasing then decreasing
haz_Ipi <- muhaz(OS.Ipi$Time, OS.Ipi$Event, bw.smooth = 3)
plot(haz_Ipi, xlab='Time (months)', main='Ipilimumab - Smoothed Hazard')

# 2.1. Standard parametric model
m_Ipi_param <- fit.models(formula = formula, data = OS.Ipi, distr = mods)

# AIC, BIC
aic_bic_summary_Ipi_param <- data.frame(AIC = sapply(m_Ipi_param$models, AIC), BIC = sapply(m_Ipi_param$models, BIC))
print(aic_bic_summary_Ipi_param)

# Visual inspection of hazard
models <- list(
  'Generalised Gamma' = m_Ipi_param$models$`Gen. Gamma`,
  'Log-Logistic' = m_Ipi_param$models$`log-Logistic`,
  'Log-Normal' = m_Ipi_param$models$`log-Normal`
)

plot_hazard(models_to_plot = models, haz = haz_Ipi, title = 'Ipilimumab - Standard Parametric Models')

# Visual inspection of survival
km_Ipi <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Ipi)
plot_survival(models_to_plot = models, km = km_Ipi, time_points = seq(0,24,0.1), title = 'Ipilimumab - Standard Parametric Models')

# 2.2 1-knot spline model
m_Ipi_spline_1_hazard <- flexsurvspline(formula = formula, data = OS.Ipi, scale = 'hazard', k = 1)
m_Ipi_spline_1_odds <- flexsurvspline(formula = formula, data = OS.Ipi, scale = 'odds', k = 1)
m_Ipi_spline_1_normal <- flexsurvspline(formula = formula, data = OS.Ipi, scale = 'normal', k = 1)

# AIC, BIC
models <- list(
  '1-knot spline hazard' = m_Ipi_spline_1_hazard,
  '1-knot spline odds' = m_Ipi_spline_1_odds,
  '1-knot spline normal' = m_Ipi_spline_1_normal
)

aic_bic_summary_Ipi_spline_1 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(Model = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Ipi_spline_1)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Ipi, title = 'Ipilimumab - 1-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Ipi, time_points = seq(0,24,0.1), title = 'Ipilimumab - 1-Knot Cubic Spline Models')

# 3. External model for both arms
haz_Scha <- muhaz(OS.Scha$Time, OS.Scha$Event, bw.smooth = 6)
plot(haz_Scha, xlab='Time (months)', main='Schadendorf Study (2015) -  Smoothed Hazard')

# 3.1. Rebased standard parametric model
rebased_time <- 13.85 # rebased_time = median follow-up of KEYNOTE006 IA2
OS.Scha.rebased <- subset(OS.Scha, Time > rebased_time)
OS.Scha.rebased$Time <-  OS.Scha.rebased$Time - rebased_time
m_Scha_rebased = fit.models(formula = formula, data = OS.Scha.rebased, distr = mods)

# AIC, BIC
aic_bic_summary_Scha_rebased <- data.frame(AIC = sapply(m_Scha_rebased$models, AIC), BIC = sapply(m_Scha_rebased$models, BIC))
print(aic_bic_summary_Scha_rebased)

# Visual inspection of hazard
models <- list(
  'Rebased Weibull' = m_Scha_rebased$models$`Weibull (AFT)`,
  'Rebased Gompertz' = m_Scha_rebased$models$Gompertz,
  'Rebased Gen Gamma' = m_Scha_rebased$models$`Gen. Gamma`,
  'Rebased Loglogistic' = m_Scha_rebased$models$`log-Logistic`,
  'Rebased Lognormal' = m_Scha_rebased$models$`log-Normal`
)

# Visual inspection of survival
plot_rebased_hazard(models_to_plot = models, haz = haz_Scha, rebased_time = rebased_time, title = 
                    'Schadendorf (2015) - Hazard of Rebased Parametric Models')

km_Scha <- survfit(Surv(Time, Event) ~ 1, data = OS.Scha)
plot_rebased_survival(models_to_plot = models, km = km_Scha, rebased_time = rebased_time, time_points = seq(0,84,by=0.1), 
                      title = 'Schadendorf (2015) - Survival of Rebased Parametric Models')

# 4. Blended method implementation
# Selected models
m_Pem_selected <- m_Pem_spline_3_normal # Internal (Pembrolizumab): 3-knot spline normal
m_Ipi_selected <- m_Ipi_param$models$`Gen. Gamma` # Internal (Ipilimumab): generalised Gamma
m_Scha_selected <- m_Scha_rebased$models$Gompertz # External (both arms): rebased Gompterz
rebased_time <- 13.85

# Parameters for blended method
t1 <- 24
t2 <- 60
a <- 5
b <- 5

# Time points
time_horizon <- 84 
dt = 0.1
t_seq <- seq(0, time_horizon, dt)

# Compute blended/internal/external hazard
h_Pem_selected <- extract_survival_hazard(m_Pem_selected, t_seq)$h
h_Ipi_selected <- extract_survival_hazard(m_Ipi_selected, t_seq)$h
h_Scha_selected <- extract_survival_hazard(m_Scha_selected, t_seq, rebased_time)$h
h_Pem_blended <- compute_blended_hazard(h_Pem_selected, h_Scha_selected, rebased_time, t1, t2, a, b, t_seq) 
h_Ipi_blended <- compute_blended_hazard(h_Ipi_selected, h_Scha_selected, rebased_time, t1, t2, a, b, t_seq)
h_Ipi_blended$est[1] <- 0

# Plot blended hazard (vs fitted internal & external hazard)
hazplot_Pem_blended <- plot_blended_hazard(h_Pem_blended, h_Pem_selected, h_Scha_selected, rebased_time, t1, t2, 'Pembrolizumab Arm')
hazplot_Ipi_blended <- plot_blended_hazard(h_Ipi_blended, h_Ipi_selected, h_Scha_selected, rebased_time, t1, t2, 'Ipilimumab Arm')
print(ggarrange(hazplot_Pem_blended, hazplot_Ipi_blended))

# 5. Comparison of extrapolated hazard and survival curves (blended method vs updated 7-year data vs TA366 base case)
# 5.1. Smoothed hazard of updated 7-year data
haz_Pem_7y <- muhaz(OS.Pem.7y$Time, OS.Pem.7y$Event, bw.smooth=6, max.time=84, n.est.grid=841)
haz_Ipi_7y <- muhaz(OS.Ipi.7y$Time, OS.Ipi.7y$Event, bw.smooth=6, max.time=84, n.est.grid=841)

# 5.2. Reproduce TA366 base case (the best standard parametric model rebased at 12 months - rebased Gompertz)
rebased_time_TA366 <- 12
OS.Scha.rebased_TA366 <- subset(OS.Scha, Time > rebased_time_TA366)
OS.Scha.rebased_TA366$Time <-  OS.Scha.rebased_TA366$Time - rebased_time_TA366
m_Scha_rebased_TA366 = fit.models(formula = formula, data = OS.Scha.rebased_TA366, distr = mods) 

aic_bic_summary_Scha_rebased_TA366 <- data.frame(AIC = sapply(m_Scha_rebased_TA366$models, AIC), BIC = sapply(m_Scha_rebased_TA366$models, BIC))
print(aic_bic_summary_Scha_rebased_TA366)

models <- list(
  'Rebased Weibull' = m_Scha_rebased_TA366$models$`Weibull (AFT)`,
  'Rebased Gompertz' = m_Scha_rebased_TA366$models$Gompertz,
  'Rebased Gen Gamma' = m_Scha_rebased_TA366$models$`Gen. Gamma`,
  'Rebased Loglogistic' = m_Scha_rebased_TA366$models$`log-Logistic`,
  'Rebased Lognormal' = m_Scha_rebased_TA366$models$`log-Normal`
)

plot_rebased_hazard(models_to_plot = models, haz = haz_Scha, rebased_time = rebased_time_TA366, title = 
                      'Schadendorf (2015) - Hazard of Parametric Models Rebased at 12 Months')

plot_rebased_survival(models_to_plot = models, km = km_Scha, rebased_time = rebased_time_TA366, time_points = seq(0,84,by=0.1), 
                      title = 'Schadendorf (2015) - Survival of Parametric Models Rebased at 12 Months')

# Extract hazard of TA366
h_TA366 <- extract_survival_hazard(m_Scha_rebased_TA366$models$Gompertz, t_seq, rebased_time_TA366)$h

# 5.3. Comparison at hazard scale (blended vs updated vs TA366)
hazplot_Pem_comparison <- plot_hazard_comparison(h_Pem_blended, haz_Pem_7y, h_TA366, rebased_time_TA366, t1, t2, 'Pembrolizumab Arm')
hazplot_Ipi_comparison <- plot_hazard_comparison(h_Ipi_blended, haz_Ipi_7y, h_TA366, rebased_time_TA366, t1, t2, 'Ipilimumab Arm')
print(ggarrange(hazplot_Pem_comparison, hazplot_Ipi_comparison))

# 5.4. Comparison at survival scale (blended vs updated vs TA366)
# KM of 7-year follow-up data
km_Pem_7y <- survfit(Surv(Time, Event)~Treatment, data=OS.Pem.7y)
km_Ipi_7y <- survfit(Surv(Time, Event)~Treatment, data=OS.Ipi.7y)

# Compute survival of blended method
S_Pem_blended <- exp(-cumsum(h_Pem_blended$est) * dt)
S_Ipi_blended <- exp(-cumsum(h_Ipi_blended$est) * dt)

# Compute survival of TA366 base case
S_Pem_rebased_time <- surv.km(OS.Pem$Time, OS.Pem$Event, rebased_time_TA366)$S.estimate
S_Ipi_rebased_time <- surv.km(OS.Ipi$Time, OS.Ipi$Event, rebased_time_TA366)$S.estimate

S_Pem_TA366_fit <- S_Pem_rebased_time * exp(-cumsum(subset(h_TA366, time >= rebased_time_TA366)$est) * dt)
S_Ipi_TA366_fit <- S_Ipi_rebased_time * exp(-cumsum(subset(h_TA366, time >= rebased_time_TA366)$est) * dt)

t_seq_km <- t_seq[t_seq < 12]
S_Pem_TA366_km <- surv.km(OS.Pem$Time, OS.Pem$Event, t_seq_km)$S.estimate
S_Pem_TA366_km[is.na(S_Pem_TA366_km)] <- 1
S_Ipi_TA366_km <- surv.km(OS.Ipi$Time, OS.Ipi$Event, t_seq_km)$S.estimate
S_Ipi_TA366_km[is.na(S_Ipi_TA366_km)] <- 1

S_Pem_TA366 <- c(S_Pem_TA366_km, S_Pem_TA366_fit)
S_Ipi_TA366 <- c(S_Ipi_TA366_km, S_Ipi_TA366_fit)

# Comparison of survival
survplot_comparison <- plot_survival_comparison(S_Pem_blended, S_Ipi_blended, km_Pem_7y, km_Ipi_7y, S_Pem_TA366, S_Ipi_TA366, rebased_time_TA366 ,t_seq)

# 6. Comparison of 7-year RMST and incremental RMST (blended vs updated vs TA366)
# Blended model
RMST_Pem_blended <- round(trapz(t_seq, S_Pem_blended) / 12, 2)
RMST_Ipi_blended <- round(trapz(t_seq, S_Ipi_blended) / 12, 2)
inc_RMST_blended <- RMST_Pem_blended - RMST_Ipi_blended

# Updated 7-year data
S_Pem_7y <- surv.km(OS.Pem.7y$Time, OS.Pem.7y$Event, t_seq)$S.estimate
S_Ipi_7y <- surv.km(OS.Ipi.7y$Time, OS.Ipi.7y$Event, t_seq)$S.estimate
S_Pem_7y[is.na(S_Pem_7y)] <- 1
S_Ipi_7y[is.na(S_Ipi_7y)] <- 1

RMST_Pem_7y <- round(trapz(t_seq, S_Pem_7y) / 12, 2)
RMST_Ipi_7y <- round(trapz(t_seq, S_Ipi_7y) / 12, 2)
inc_RMST_7y <- RMST_Pem_7y - RMST_Ipi_7y

# TA366 base case (0-12 mo Kaplan-Meier + 12-84 mo model fit)
RMST_Pem_TA366 <- round(trapz(t_seq, S_Pem_TA366) / 12, 2)
RMST_Ipi_TA366 <- round(trapz(t_seq, S_Ipi_TA366) / 12, 2)
inc_RMST_TA366 <- RMST_Pem_TA366 - RMST_Ipi_TA366

# RMST table
RMST_table <- matrix(
  c(RMST_Pem_blended, RMST_Ipi_blended, inc_RMST_blended,
    RMST_Pem_7y, RMST_Ipi_7y, inc_RMST_7y,
    RMST_Pem_TA366, RMST_Ipi_TA366, inc_RMST_TA366),
  nrow = 3, byrow = TRUE
)

colnames(RMST_table) <- c('Pembrolizumab RMST 7y', 'Ipilimumab RMST 7y', 'Incremental RMST 7y')
rownames(RMST_table) <- c('Blended Method', '7-Year Updated', 'TA366 Base Case')

write.csv(RMST_table, "tables/RMST_base_case.csv", row.names = T)
