# Function to plot hazard 
plot_hazard <- function(models_to_plot, haz, title){
  haz_df <- data.frame(Time = haz$est.grid, Hazard = haz$haz.est)
  hazard_plots <- lapply(names(models_to_plot), function(model_name) {
    haz_data <- as.data.frame(summary(models_to_plot[[model_name]], type = "hazard"))
    ggplot() +
      geom_line(data = haz_df, aes(x = Time, y = Hazard, color = "Smoothed Hazard"), linewidth = 1) +
      geom_line(data = haz_data, aes(x = time, y = est, color = "Model Fit"), linewidth = 1) +
      geom_ribbon(data = haz_data, aes(x = time, y = est, ymin = lcl, ymax = ucl), 
                  alpha = 0.1, fill = "#377EB8") +
      scale_colour_manual(name = "", values = c("Smoothed Hazard" = "black", "Model Fit" = "#377EB8")) +
      geom_vline(xintercept = 13.85, linetype = "dashed") +
      xlab("Time (months)") + ylab("Hazard") + ylim(0, 0.085) + theme_classic() +
      theme(legend.position = "top", legend.direction = "horizontal") +
      ggtitle(model_name)
  })
  plot <- annotate_figure(ggarrange(plotlist = hazard_plots, nrow = 1), top = text_grob(title, face = "bold"))
  print(plot)
}

# Function to plot survival
plot_survival <- function(models_to_plot, km, time_points, title){
  km_df <- data.frame(time = km$time, est = km$surv, Model = "Kaplan-Meier")
  survival_estimates <- lapply(names(models_to_plot), function(model_name) {
    surv_data <- as.data.frame(summary(models_to_plot[[model_name]], type = "survival", t = time_points))[,c('time','est')]
    surv_data$Model <- model_name  
    return(surv_data)
  })
  surv_df <- rbind(do.call(rbind, survival_estimates), km_df)
  brewer_palette <- brewer.pal(length(models_to_plot), "Set1")[1:length(models_to_plot)]
  model_colors <- setNames(brewer_palette, names(models_to_plot))
  model_colors <- c("Kaplan-Meier" = "black", model_colors)
  plot <- ggplot(surv_df, aes(x = time, y = est, color = Model)) +
    geom_line(linewidth = 1) +  
    scale_color_manual(values = model_colors) +
    xlab("Time (months)") + ylab("Overall Survival") + 
    ylim(0,1) + scale_x_continuous(limits = c(0,24), breaks = seq(0,24,by=4)) +
    ggtitle(title) + theme_classic() + theme(legend.position = "top")
  print(plot)
}

# Function to plot rebased hazard
plot_rebased_hazard <- function(models_to_plot, haz, rebased_time, title) {
  haz_df <- data.frame(Time = haz$est.grid, Hazard = haz$haz.est)
  hazard_plots <- lapply(names(models_to_plot), function(model_name) {
    haz_data <- as.data.frame(summary(models_to_plot[[model_name]], type = "hazard"))
    haz_data$time <- haz_data$time + rebased_time
    
    plot <- ggplot() +
      geom_line(data = haz_df, aes(x = Time, y = Hazard, color = "Smoothed Hazard"), linewidth = 1) +
      geom_line(data = haz_data, aes(x = time, y = est, color = "Model Fit"), linewidth = 1) +
      geom_ribbon(data = haz_data, aes(x = time, y = est, ymin = lcl, ymax = ucl), 
                  alpha = 0.1, fill = "#377EB8") +
      scale_colour_manual(name = "", values = c("Smoothed Hazard" = "black", "Model Fit" = "#377EB8")) +
      geom_vline(xintercept = rebased_time, linetype = "dashed") +
      xlab("Time (months)") + ylab("Hazard") + ylim(0, 0.09) + theme_classic() +
      theme(legend.position = "top", legend.direction = "horizontal") +
      ggtitle(model_name)
  })
  num_models <- length(models_to_plot)
  if (num_models <= 3) {
    nrow <- 1
    ncol <- num_models  
  } else {
    nrow <- 2
    ncol <- ceiling(num_models / 2) 
  }
  plot <- annotate_figure(ggarrange(plotlist = hazard_plots, nrow = nrow, ncol = ncol), 
                          top = text_grob(title, face = "bold"))
  print(plot)
}

# Function to plot rebased survival
plot_rebased_survival <- function(models_to_plot, km, rebased_time, time_points, title){
  km_df <- data.frame(time = km$time, est = km$surv, Model = "Kaplan-Meier")
  S_rebased_time <- surv.km(OS.Scha$Time, OS.Scha$Event, rebased_time)$S.estimate
  survival_estimates <- lapply(names(models_to_plot), function(model_name) {
    surv_data <- as.data.frame(summary(models_to_plot[[model_name]], type = "survival", t = time_points - rebased_time))[,c('time','est')]
    surv_data$time <- surv_data$time + rebased_time  
    surv_data$est <- S_rebased_time * surv_data$est
    surv_data$Model <- model_name
    surv_data <- surv_data[surv_data$time >= rebased_time, ] # Keep only values where time >= rebased_time
    return(surv_data)
  })
  surv_df <- rbind(do.call(rbind, survival_estimates), km_df)
  
  num_models <- length(models_to_plot)
  brewer_palette <- brewer.pal(max(3, num_models), "Set1")[1:num_models]
  model_colors <- setNames(brewer_palette, names(models_to_plot))
  model_colors <- c("Kaplan-Meier" = "black", model_colors)
  
  plot <- ggplot(surv_df, aes(x = time, y = est, color = Model)) +
    geom_line(linewidth = 1) +  
    geom_vline(xintercept = rebased_time, linetype = "dashed") +
    scale_color_manual(values = model_colors) +
    xlab("Time (months)") + ylab("Overall Survival") + 
    ylim(0,1) + scale_x_continuous(limits = c(0,84), breaks = seq(0,84,by=4)) +
    ggtitle(title) + theme_classic() + theme(legend.position = "top")
  
  print(plot)
}

# Function to extract survival and hazard of selected models
extract_survival_hazard <- function(model, t_seq, rebased_time = NULL) {
  if (!is.null(rebased_time)){
    S_data <- as.data.frame(summary(model, t = t_seq - rebased_time))
    h_data <- as.data.frame(summary(model, type = "hazard", t = t_seq - rebased_time))
    
    # Adjust time
    S_data$time <- S_data$time + rebased_time
    h_data$time <- h_data$time + rebased_time
    
    # Adjust survival probabilities for rebased models
    S_rebased_time <- surv.km(OS.Scha$Time, OS.Scha$Event, rebased_time)$S.estimate
    S_data$est <- S_rebased_time * S_data$est
    
    # Keep only time points where time >= rebased_time
    S_data <- S_data[S_data$time >= rebased_time, ]
  }
  else{
    S_data <- as.data.frame(summary(model, t = t_seq))
    h_data <- as.data.frame(summary(model, type = "hazard", t = t_seq))
  }
  
  return(list(S = S_data, h = h_data))
}

# Function to plot survival of selected models and KMs
plot_survival_selected <- function(m_Pem, m_Ipi, m_Scha, rebased_time, t_seq){
  # Extract fitted survival
  S_Pem_selected <- extract_survival_hazard(m_Pem_selected, t_seq)$S
  S_Ipi_selected <- extract_survival_hazard(m_Ipi_selected, t_seq)$S
  S_Scha_selected <- extract_survival_hazard(m_Scha_selected, t_seq, rebased_time)$S
  
  # Combine fitted survival
  model_survival_selected_df <- rbind(
    data.frame(time = S_Pem_selected$time, est = S_Pem_selected$est, Group = "Pembrolizumab", Type = "Selected Model Fit"),
    data.frame(time = S_Ipi_selected$time, est = S_Ipi_selected$est, Group = "Ipilimumab", Type = "Selected Model Fit"),
    data.frame(time = S_Scha_selected$time, est = S_Scha_selected$est, Group = "Schadendorf", Type = "Selected Model Fit")
  )
  
  # Combine Kaplan-Meier data
  OS.int.ext <- list("Pembrolizumab" = OS.Pem, "Ipilimumab" = OS.Ipi, "Schadendorf" = OS.Scha)
  km_int_ext <- lapply(OS.int.ext, function(data) survfit(Surv(Time, Event) ~ Treatment, data = data))
  km_int_ext_df <- do.call(rbind, lapply(names(km_int_ext), function(name) {
    km_fit <- km_int_ext[[name]]
    data.frame(time = km_fit$time, surv = km_fit$surv, Group = name, Type = "Kaplan-Meier")
  }))
  
  # Plot 
  ggplot() +
    geom_line(data = km_int_ext_df, aes(x = time, y = surv, linetype = Type, color = Group), linewidth = 1) +
    geom_line(data = model_survival_selected_df, aes(x = time, y = est, linetype = Type, color = Group), linewidth = 0.8) +
    scale_linetype_manual(name = "Profile", values = c("Kaplan-Meier" = "solid", "Selected Model Fit" = "dashed")) +
    scale_color_brewer(palette = "Set1", name = "Data") +
    scale_x_continuous(name = "Time (months)", breaks = seq(0, 84, 12), limits = c(0, 84)) +
    scale_y_continuous(name = "Overall Survival", limits = c(0, 1)) +
    theme_classic() +
    ggtitle('Fitted Survival of Selected Models')
}

# Function to compute blended hazard 
compute_blended_hazard <- function(h_int, h_ext, rebased_time, t1, t2, a, b, t_seq) {
  weight <- pbeta((t_seq - t1) / (t2 - t1), shape1 = a, shape2 = b)
  h_blended_est <- (1 - weight) * h_int$est + weight * h_ext$est
  
  return(data.frame(time = t_seq, est = h_blended_est))
}

# Function to plot blended hazard (vs fitted internal/external hazard)
plot_blended_hazard <- function(h_blended, h_int, h_ext, rebased_time, t1, t2, title) {
  ggplot() +
    geom_line(data = h_blended, aes(x = time, y = est, color = "Blended Hazard"), linewidth = 1) + 
    geom_line(data = h_int, aes(x = time, y = est, color = "Fitted Internal Hazard"), linewidth = 1, linetype = 'dashed') +
    geom_line(data = subset(h_ext, time >= rebased_time), aes(x = time, y = est, color = "Fitted External Hazard"), linewidth = 1, linetype = 'dashed') +
    geom_vline(xintercept = c(t1, t2), linetype="dashed", color = "grey") +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) + ylim(0,0.065) +
    scale_color_brewer(palette = "Set1", name = "Model") +
    theme_classic() + theme(legend.position = c(0.75, 0.8), text = element_text(size = 12), legend.key.width = unit(1, "cm")) +
    labs(title = title, x = "Time (months)", y = "Hazard")
}

# Function to compare blended hazard vs smoothed updated hazard vs TA366 estimated hazard
plot_hazard_comparsion <- function(h_blended, haz_updated, h_TA366, rebased_time, t1, t2, title) {
  ggplot() +
    geom_line(data = h_blended, aes(x = time, y = est, color = "Blended Model"), linewidth = 1) + 
    geom_line(aes(haz_updated$est.grid, haz_updated$haz.est, color="Updated Data"), linewidth = 1) +
    geom_line(data = subset(h_TA366, time >= rebased_time_TA366), aes(x = time, y = est, color = "TA366 Base Case"), linewidth = 1) +
    geom_vline(xintercept = c(t1, t2), color = "grey", linetype = 'dashed') +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) + ylim(0, 0.065) +
    scale_color_brewer(palette = "Set1", name = "Model") +
    theme_classic() + theme(legend.position = c(0.75, 0.8), text = element_text(size = 12), legend.key.width = unit(1, "cm")) +
    labs(title = title, x = "Time (months)", y = "Hazard")
}

# Function to compare blended survival vs updated KM vs TA366 estimated survival
plot_survival_comparsion <- function(S_Pem_blended, S_Ipi_blended, km_Pem_updated, km_Ipi_updated, S_Pem_TA366, S_Ipi_TA366, rebased_time_TA366, t_seq){
  plot <- ggplot() +
    # Blended survival estimates
    geom_line(aes(x = t_seq, y = S_Pem_blended, color = "Blended Model", linetype = "Pembrolizumab"), linewidth = 1) +
    geom_line(aes(x = t_seq, y = S_Ipi_blended, color = "Blended Model", linetype = "Ipilimumab"), linewidth = 1) +
    # Updated Kaplan-Meier
    geom_line(aes(x = km_Pem_7y$time, y = km_Pem_7y$surv, color = 'Updated Kaplan-Meier', linetype = "Pembrolizumab"), linewidth = 1) +
    geom_line(aes(x = km_Ipi_7y$time, y = km_Ipi_7y$surv, color = 'Updated Kaplan-Meier', linetype = "Ipilimumab"), linewidth = 1) +
    # TA366 base case 
    geom_line(aes(x = t_seq[t_seq >= rebased_time_TA366], y = S_Pem_TA366, color = 'TA366 Base Case', linetype = "Pembrolizumab"), linewidth = 1) +
    geom_line(aes(x = t_seq[t_seq >= rebased_time_TA366], y = S_Ipi_TA366, color = 'TA366 Base Case', linetype = "Ipilimumab"), linewidth = 1) +
    scale_x_continuous(name = "Time (months)", breaks = seq(0, 84, 12), limits = c(0, 84)) +
    scale_y_continuous(name = "Overall Survival", limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1", name = "Model") +
    scale_linetype_manual(name = "Treatment", values = c("Pembrolizumab" = "solid", "Ipilimumab" = "dashed")) +
    theme_classic() +
    theme(legend.direction = "vertical", legend.box = "vertical", legend.position = c(0.8, 0.8), legend.key.width = unit(1, "cm"))
  
  print(plot)
}