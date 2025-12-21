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
formula <- Surv(Time, Event) ~ 1
mods <- c('exp', 'weibull', 'gompertz', 'gengamma', 'loglogistic', 'lognormal')
m_Pem_param <- fit.models(formula = formula, data = OS.Pem, distr = mods)

models <- list(
  'Generalised Gamma' = m_Pem_param$models$`Gen. Gamma`,
  'Log-Logistic' = m_Pem_param$models$`log-Logistic`,
  'Log-Normal' = m_Pem_param$models$`log-Normal`
)

# Visual inspection of hazard plot
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - Standard Parametric Models')

# Visual inspection of survival plot
km_Pem <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Pem)
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,20,0.1), title = 'Pembrolizumab - Standard Parametric Models')

# AIC, BIC
aic_bic_summary_Pem_param <- data.frame(AIC = sapply(m_Pem_param$models, AIC), BIC = sapply(m_Pem_param$models, BIC))

print(aic_bic_summary_Pem_param)

# 1.2. 1-knot cubic spline model
m_Pem_spline_1_hazard <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'hazard', k = 1)
m_Pem_spline_1_odds <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'odds', k = 1)
# Error in optim using 'BFGS', thus changing to 'Nelder-Mead'
m_Pem_spline_1_normal <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'normal', k = 1, method = "Nelder-Mead")

models <- list(
  '1-knot spline hazard' = m_Pem_spline_1_hazard,
  '1-knot spline odds' = m_Pem_spline_1_odds,
  '1-knot spline normal' = m_Pem_spline_1_normal
)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - 1-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,20,0.1), title = 'Pembrolizumab - 1-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Pem_spline_1 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_1)

# 1.3. 2-knot cubic spline model
m_Pem_spline_2_hazard <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'hazard', k = 2)
m_Pem_spline_2_odds <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'odds', k = 2)
m_Pem_spline_2_normal <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'normal', k = 2)

models <- list(
  '2-knot spline hazard' = m_Pem_spline_2_hazard,
  '2-knot spline odds' = m_Pem_spline_2_odds,
  '2-knot spline normal' = m_Pem_spline_2_normal
)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - 2-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0,20,0.1), title = 'Pembrolizumab - 2-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Pem_spline_2 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_2)

# 1.4. 3-knot cubic spline model
m_Pem_spline_3_hazard <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'hazard', k = 3)
m_Pem_spline_3_odds <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'odds', k = 3)
m_Pem_spline_3_normal <- flexsurvspline(formula = formula, data = OS.Pem, scale = 'normal', k = 3)

models <- list(
  '3-knot spline hazard' = m_Pem_spline_3_hazard,
  '3-knot spline odds' = m_Pem_spline_3_odds,
  '3-knot spline normal' = m_Pem_spline_3_normal
)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Pem, title = 'Pembrolizumab - 3-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Pem, time_points = seq(0, 20, 0.1), title = 'Pembrolizumab - 3-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Pem_spline_3 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Pem_spline_3)
# write.csv(rbind(aic_bic_summary_Pem_param, aic_bic_summary_Pem_spline_1, aic_bic_summary_Pem_spline_2, aic_bic_summary_Pem_spline_3), 'tables/model selection/aic_bic_int_pem_selection.csv', row.names = T)

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
plot_hazard(models_to_plot = models, haz = haz_Ipi, title = 'Ipilimumab - Standard Parametric Models')

# Visual inspection of survival
km_Ipi <- survfit(Surv(Time, Event) ~ Treatment, data = OS.Ipi)
plot_survival(models_to_plot = models, km = km_Ipi, time_points = seq(0,20,0.1), title = 'Ipilimumab - Standard Parametric Models')

# AIC, BIC
aic_bic_summary_Ipi_param <- data.frame(AIC = sapply(m_Ipi_param$models, AIC), BIC = sapply(m_Ipi_param$models, BIC))

print(aic_bic_summary_Ipi_param)

# 2.2 1-knot spline model
m_Ipi_spline_1_hazard <- flexsurvspline(formula = formula, data = OS.Ipi, scale = 'hazard', k = 1)
m_Ipi_spline_1_odds <- flexsurvspline(formula = formula, data = OS.Ipi, scale = 'odds', k = 1)
m_Ipi_spline_1_normal <- flexsurvspline(formula = formula, data = OS.Ipi, scale = 'normal', k = 1)

models <- list(
  '1-knot spline hazard' = m_Ipi_spline_1_hazard,
  '1-knot spline odds' = m_Ipi_spline_1_odds,
  '1-knot spline normal' = m_Ipi_spline_1_normal
)

# Visual inspection of hazard
plot_hazard(models_to_plot = models, haz = haz_Ipi, title = 'Ipilimumab - 1-Knot Cubic Spline Models')

# Visual inspection of survival
plot_survival(models_to_plot = models, km = km_Ipi, time_points = seq(0, 20, 0.1), title = 'Ipilimumab - 1-Knot Cubic Spline Models')

# AIC, BIC
aic_bic_summary_Ipi_spline_1 <- do.call(rbind, lapply(names(models), function(model_name) {
  data.frame(row.names = model_name, AIC = AIC(models[[model_name]]), BIC = BIC(models[[model_name]]))
}))

print(aic_bic_summary_Ipi_spline_1)
# write.csv(rbind(aic_bic_summary_Ipi_param, aic_bic_summary_Ipi_spline_1), 'tables/model selection/aic_bic_int_ipi_selection.csv', row.names = T)

# 3. External model for both arms
haz_Scha <- muhaz(OS.Scha$Time, OS.Scha$Event, bw.smooth = 6)
plot(haz_Scha, xlab='Time (months)', main='Schadendorf Study (2015) -  Smoothed Hazard')
abline(v = 13.85, col = 'grey', lty = 2)

# 3.1. Rebased standard parametric model
rebased_time <- 13.85 # rebased_time = median follow-up of KEYNOTE006 IA2
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
plot_rebased_hazard(models_to_plot = models, haz = haz_Scha, rebased_time = rebased_time, title = 
                    'Schadendorf (2015) - Hazard of Parametric Models Rebased at 13.84 Months')

# Visual inspection of survival
km_Scha <- survfit(Surv(Time, Event) ~ 1, data = OS.Scha)
plot_rebased_survival(models_to_plot = models, km = km_Scha, rebased_time = rebased_time, time_points = seq(0, 84, by = 0.1), 
                      title = 'Schadendorf (2015) - Survival of Parametric Models Rebased at 13.84 Months')

# AIC, BIC
aic_bic_summary_Scha_rebased <- data.frame(AIC = sapply(m_Scha_rebased$models, AIC), BIC = sapply(m_Scha_rebased$models, BIC))
print(aic_bic_summary_Scha_rebased)
# write.csv(aic_bic_summary_Scha_rebased, 'tables/model selection/aic_bic_ext_selection.csv', row.names = T)

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
h_Pem_blended <- compute_blended_hazard(h_Pem_selected, h_Scha_selected, t1, t2, a, b, t_seq) 
h_Ipi_blended <- compute_blended_hazard(h_Ipi_selected, h_Scha_selected, t1, t2, a, b, t_seq)
h_Ipi_blended$est[1] <- 0

# Plot blended hazard (vs fitted internal & external hazard)
hazplot_Pem_blended <- plot_blended_hazard(h_Pem_blended, h_Pem_selected, h_Scha_selected, rebased_time, t1, t2, 'Pembrolizumab')
hazplot_Ipi_blended <- plot_blended_hazard(h_Ipi_blended, h_Ipi_selected, h_Scha_selected, rebased_time, t1, t2, 'Ipilimumab')
hazplot_blended <- ggarrange(hazplot_Pem_blended, hazplot_Ipi_blended)

hazplot_blended

# 5. Comparison of extrapolated hazard and survival curves (blended method vs updated 7-year data vs TA366 base case)
# 5.1. Reproduce TA366 base case (the best standard parametric model rebased at 12 months - rebased Gompertz)
rebased_time_TA366 <- 12
OS.Scha.rebased_TA366 <- subset(OS.Scha, Time > rebased_time_TA366)
OS.Scha.rebased_TA366$Time <-  OS.Scha.rebased_TA366$Time - rebased_time_TA366
m_Scha_rebased_TA366 = fit.models(formula = formula, data = OS.Scha.rebased_TA366, distr = mods) 

models <- list(
  'Rebased Weibull' = m_Scha_rebased_TA366$models$`Weibull (AFT)`,
  'Rebased Gompertz' = m_Scha_rebased_TA366$models$Gompertz,
  'Rebased Gen Gamma' = m_Scha_rebased_TA366$models$`Gen. Gamma`,
  'Rebased Loglogistic' = m_Scha_rebased_TA366$models$`log-Logistic`,
  'Rebased Lognormal' = m_Scha_rebased_TA366$models$`log-Normal`
)

hazplot_TA366_selection <- plot_rebased_hazard(models_to_plot = models, haz = haz_Scha, rebased_time = rebased_time_TA366, title = 
                      'Schadendorf (2015) - Hazard of Parametric Models Rebased at 12 Months')

survplot_TA366_selection <- plot_rebased_survival(models_to_plot = models, km = km_Scha, rebased_time = rebased_time_TA366, time_points = seq(0,84,by=0.1), 
                      title = 'Schadendorf (2015) - Survival of Parametric Models Rebased at 12 Months')

aic_bic_summary_Scha_rebased_TA366 <- data.frame(AIC = sapply(m_Scha_rebased_TA366$models, AIC), BIC = sapply(m_Scha_rebased_TA366$models, BIC))
# write.csv(aic_bic_summary_Scha_rebased_TA366, 'tables/model selection/aic_bic_TA366_selection.csv', row.names = T)

# 5.2. Comparison at hazard scale (blended vs updated vs TA366)
haz_Pem_7y <- muhaz(OS.Pem.7y$Time, OS.Pem.7y$Event, bw.smooth=6, max.time=84, n.est.grid=841)
haz_Ipi_7y <- muhaz(OS.Ipi.7y$Time, OS.Ipi.7y$Event, bw.smooth=6, max.time=84, n.est.grid=841)
h_TA366 <- extract_survival_hazard(m_Scha_rebased_TA366$models$Gompertz, t_seq, rebased_time_TA366)$h

hazplot_Pem_comparison <- plot_hazard_comparison(h_Pem_blended, haz_Pem_7y, h_TA366, rebased_time_TA366, t1, t2, 'Pembrolizumab')
hazplot_Ipi_comparison <- plot_hazard_comparison(h_Ipi_blended, haz_Ipi_7y, h_TA366, rebased_time_TA366, t1, t2, 'Ipilimumab')
hazplot_comparison <-  ggarrange(hazplot_Pem_comparison, hazplot_Ipi_comparison)

hazplot_comparison

# 5.3. Comparison at survival scale (blended vs updated vs TA366)
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
print(RMST_table)

# 7. Comparison of HR (blended vs updated vs TA366)
HR_blend_df <- merge(h_Pem_blended, h_Ipi_blended, by = "time", suffixes = c("_Pem", "_Ipi"))
HR_blend_df$HR <- HR_blend_df$est_Pem / HR_blend_df$est_Ipi
HR_blend_df_sub <- subset(HR_blend_df, time > 0.2)

haz_Pem_7y <- muhaz(OS.Pem.7y$Time, OS.Pem.7y$Event, bw.smooth=24, max.time=84, n.est.grid=841)
haz_Ipi_7y <- muhaz(OS.Ipi.7y$Time, OS.Ipi.7y$Event, bw.smooth=24, max.time=84, n.est.grid=841)

# Define function for one bootstrap iteration
bootstrap_HR <- function(data_pem, data_ipi, time_grid) {
  # Resample with replacement
  resample_pem <- data_pem[sample(nrow(data_pem), replace = TRUE), ]
  resample_ipi <- data_ipi[sample(nrow(data_ipi), replace = TRUE), ]
  
  # Recalculate smoothed hazard using muhaz
  haz_pem <- muhaz(resample_pem$Time, resample_pem$Event, bw.smooth = 24, 
                   max.time = 84, n.est.grid = length(time_grid))
  haz_ipi <- muhaz(resample_ipi$Time, resample_ipi$Event, bw.smooth = 24, 
                   max.time = 84, n.est.grid = length(time_grid))
  
  # Interpolate hazard estimates if necessary
  h_pem <- approx(haz_pem$est.grid, haz_pem$haz.est, xout = time_grid, rule = 1)$y
  h_ipi <- approx(haz_ipi$est.grid, haz_ipi$haz.est, xout = time_grid, rule = 1)$y
  
  # Calculate HR
  HR <- h_pem / h_ipi
  return(HR)
}

set.seed(123)  # for reproducibility

B <- 500  # Number of bootstrap iterations
time_grid <- seq(0, 84, length.out = 841)  # Same as original grid

# Run bootstrap
boot_HR_mat <- replicate(B, bootstrap_HR(OS.Pem.7y, OS.Ipi.7y, time_grid))

# boot_HR_mat: rows = time points, columns = bootstrap replicates

# Bootstrap percentiles
HR_smoothed_df <- data.frame(
  time = haz_Pem_7y$est.grid,
  HR = haz_Pem_7y$haz.est / haz_Ipi_7y$haz.est,
  lower = apply(boot_HR_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  upper = apply(boot_HR_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
)

HRplot_comparison <- ggplot() + 
  geom_line(data = HR_blend_df_sub, aes(x = time, y = HR, color = 'Blended hazard method'), linewidth = 1) +
  geom_line(data = HR_smoothed_df, aes(x = time, y = HR, color = 'Smoothed HR of updated data'), linewidth = 1) +
  geom_ribbon(data = HR_smoothed_df, aes(x = time, ymin = lower, ymax = upper), fill = "#90EE90", alpha = 0.2) +
  geom_segment(aes(x = 12, xend = 84, y = 1, yend = 1, color = 'Piecewise method in TA366 base case'), linewidth = 1) +
  geom_hline(yintercept = 1, linetype='dashed', color = 'grey') +
  geom_vline(xintercept = c(t1, t2), linetype='dashed', color = 'grey') +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) + ylim(0, 2.5) +
  scale_color_brewer(palette = 'Set1', name = '') +
  theme_classic() + theme(text = element_text(size = 12), legend.position = c(0.8, 0.85), legend.key.width = unit(1, 'cm')) +
  labs(x = 'Time (months)', y = 'HR')

#Save output figures and tables
if (t1 == 24 & t2 == 60 & a == 5 & b == 5){
  ggsave(paste0('figures/base case/hazplot_blended_', t1, '_', t2, '_', a, '_', b, '.png'),
         hazplot_blended, width = 8, height = 5, units = 'in')
  ggsave(paste0('figures/base case/hazplot_comparison_', t1, '_', t2, '_', a, '_', b, '.png'),
         hazplot_comparison, width = 8, height = 5, units = 'in')
  ggsave(paste0('figures/base case/survplot_comparison_', t1, '_', t2, '_', a, '_', b, '.png'),
         survplot_comparison, width = 6, height = 5, units = 'in')
  ggsave(paste0('figures/base case/hrplot_comparison_', t1, '_', t2, '_', a, '_', b, '.png'),
         HRplot_comparison, width = 6, height = 5, units = 'in')
  write.csv(RMST_table, paste0('tables/base case/RMST_', t1, '_', t2, '_', a, '_', b, '.csv'), row.names = T)
} else{
  ggsave(paste0('figures/sensitivity analysis/hazplot_blended_', t1, '_', t2, '_', a, '_', b, '.png'),
         hazplot_blended, width = 8, height = 5, units = 'in')
  ggsave(paste0('figures/sensitivity analysis/hazplot_comparison_', t1, '_', t2, '_', a, '_', b, '.png'),
         hazplot_comparison, width = 8, height = 5, units = 'in')
  ggsave(paste0('figures/sensitivity analysis/survplot_comparison_', t1, '_', t2, '_', a, '_', b, '.png'),
         survplot_comparison, width = 6, height = 5, units = 'in')
  ggsave(paste0('figures/sensitivity analysis/hrplot_comparison_', t1, '_', t2, '_', a, '_', b, '.png'),
         HRplot_comparison, width = 6, height = 5, units = 'in')
  write.csv(RMST_table, paste0('tables/sensitivity analysis/RMST_', t1, '_', t2, '_', a, '_', b, '.csv'), row.names = T)
}
