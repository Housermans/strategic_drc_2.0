


# # Plotten van de curves met ggplot2
# ggplot(experiment_df,
#        aes(x = conc_condition,
#            y = mean_GR)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower_bound_95,
#                     ymax = upper_bound_95)) +
#   facet_wrap(~condition) +
#   labs(x = "Concentratie",
#        y = "Gemiddelde GR") +
#   stat_smooth(method=drm, fct=LL.4(), se=FALSE)
# 
# # Maak een lege lijst om de fit objecten op te slaan
# fit_list <- list()
# 
# # Loop over de unieke waarden van condition
# for (cond in unique(experiment_df$condition)) {
#   # Maak een subset van de data voor elke condition
#   sub_data <- experiment_df[experiment_df$condition == cond, ]
#   # Pas de drm functie toe op de subset data
#   fit <- drc::drm(mean_GR ~ conc_condition, data = sub_data, fct = LL.4())
#   # Sla het fit object op in de lijst met de naam van de condition
#   fit_list[[cond]] <- fit
# }
# 
# # Bekijk de lijst met fit objecten
# fit_list
# 
# # Write a function that makes a plot for a given condition
# plot_condition <- function(condition) {
#   
#   # Subset the data frame by condition using filter()
#   df <- filter(experiment_df,
#                condition == condition)
#   
#   # Make a plot with log scale x axis and adjusted range
#   p <- ggplot(df,
#               aes(x = conc_condition,
#                   y = mean_GR)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = lower_bound_95,
#                       ymax = upper_bound_95)) +
#     labs(x = "Concentration",
#          y = "Mean GR",
#          title = paste("Condition", condition)) +
#     stat_smooth(method=drm,fct=LL.4(),se=FALSE) +
#     
#     # Add log scale transformation and limits
#     scale_x_log10(limits = c(0.001,100))
#   
#   # Return the plot
#   return(p)
# }
# 
# # Loop over the unique values of condition
# for (c in unique(experiment_df$condition)) {
#   
#   # Clear the plot device
#   dev.off()
#   
#   # Call the function and assign the result to p_c
#   p_c <- plot_condition(c)
#   
#   # Print out or save the plot as desired
#   print(p_c)
#   ggsave(p_c,file.path(plot_dir, paste("p_",c,"png",sep="")))
# }