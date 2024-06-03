#---------------------------------------------#
#                   SET UP
#---------------------------------------------#
rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

# PATH
PATH <- "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/olink/"

# read in mediation results
mediation_results <- openxlsx::read.xlsx(paste0(PATH, "mediation-olink/mediation-results-CM-olink-MM-adjusted-2024-02-11.xlsx"), rowNames = T)

#---------------------------------------------#
#                     FORMAT
#---------------------------------------------#
# This script plots the MR mediation results (a path, b path and indirect effect)

# reduce opacity of grid lines (default is 255)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

#---------------------------------#
#    DEFINE PLOT (indirect effect)
#---------------------------------#
# plot indirect effects from product of coef
# define indirect effect plotting function
indirect_effect_plot <- function(data, trait) {
  # save markers as factor and arrange accordingly
  data$marker <- as.factor(data$marker) 
  data <- data %>% arrange(marker) # makes sure order of pvalues and factor names are the same (important as the order of levels in the factor will determine the order in the plot)
  
  # # add colours for sig. effects
  # ind <- which(data$ab.pval < 0.05)
  # my_colors <- rep("grey", nrow(data))
  # my_colors[ind] <- "black" # Significant effects in black
  data$marker_numeric <- data$marker
  
  # choose colours for forest plot
  my_colors <- rep("#D9CDED", nrow(data)) # can use blue instead of purple or green "#a9def9"
  # add colours for lipids
  ind <- which(data$category == "Lipids")
  my_colors[ind] <-"#EDE981" #  lipid effects in yellow
  # add colour for cortisol
  ind <- which(data$category == "Cortisol")
  my_colors[ind] <- "#d0f4de"  # "#EBCADF" # Cortisol effects in green
  # add colours for glycaemic
  ind <- which(data$category == "Glycaemic")
  my_colors[ind] <- "#ff99c8" # "#EBCADF" # Glycaemic effects in pink
  
  forest3 <- ggplot(
    data = data,
    aes(x = factor(marker_numeric), y = ab, ymin = ab.lCI, ymax = ab.hCI)
  ) +
    geom_errorbar(aes(col = marker), color = "darkgrey", width = 0.4  # Set the width of the error bars
    ) +
    geom_pointrange(aes(col = marker), size = 1.6) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    xlab(" ") + # was "Indirect effect"
    ylab("beta (95% CI)") +
    theme_minimal() +
    theme(
      panel.background = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      strip.text.y = element_text(face = "bold", size = 20),  
      panel.grid.major.y = element_line(color = col_grid, size = 0.5),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.position = "none",
      axis.text = element_text(face = "plain", size = 16),
      axis.title = element_text(face = "plain", size = 16),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ) + scale_color_manual(values = my_colors ) +
    coord_flip() + 
    # ggtitle(paste0("CM-", trait)) + # add title here
    theme(plot.title = element_text(size = 16, face = "bold")
    )  +
 scale_x_discrete(labels = c("CSF.1", "LIF.R", "C-reactive protein", expression("HbA1c"^c), expression("HbA1c"^b), expression("HDL cholesterol"^a), expression("HDL cholesterol"^b), expression("LDL cholesterol"^c), expression("Triglycerides"^a), expression("Triglycerides"^b)))
  
  return(forest3)
}

# add same names as for a and b paths
mediation_results$marker <- sub("-1", "", row.names(mediation_results))
row.names(mediation_results) <- mediation_results$marker # add mediator name as row name (remove '-1' from marker names)

# define data format function
my_data_format <- function(data){
  
  # add marker names
  data$marker <- row.names(data)
  data$ab.lCI <- data$ab -1.96*data$ab.se
  data$ab.hCI <- data$ab + 1.96*data$ab.se
  
  # order alphabetically based on marker name
  data <- data[order(data$marker),]
  return(data)
}

# format data for plotting
mediation_results_formatted <- my_data_format(data=mediation_results)

#--------------------------------------#
#   PLOT indirect effect for all traits
#--------------------------------------#
# combine indirect effects into one plot for all traits (includes 'other traits' and 'olink' markers)

# read in mediation results for other traits
mediation_results_other_traits <- openxlsx::read.xlsx(paste0(PATH, "../other-traits/mediation-other-traits/mediation-results-CM-other-traits-MM-adjusted-2024-02-10.xlsx"), rowNames = T)

# add same names as for a and b paths
mediation_results_other_traits$marker <- row.names(mediation_results_other_traits)

# format data for plotting
mediation_results_other_traits_formatted <- my_data_format(data=mediation_results_other_traits)

# combine olink with other traits
all_traits <- rbind(mediation_results_formatted, mediation_results_other_traits_formatted)

# add category variable (used for assigning colours in plot)
all_traits$category <- NA
all_traits$category <- ifelse(rownames(all_traits) %in% c("CSF.1", "LIF.R", "CRP"), "Inflammatory",
                              ifelse(rownames(all_traits) %in% c("HbA1c-ebi-a-GCST90014006", "HbA1c-ieu-b-4842"), "Glycaemic",
                                     ifelse(rownames(all_traits) %in% c("HDL-ieu-a-299", "HDL-ieu-b-109", "LDL-ebi-a-GCST90018961", "TG-ieu-a-302", "TG-ieu-b-111"), "Lipids", NA)))

# Change the levels of the marker variable to desired order in plot
all_traits$marker <- factor(all_traits$marker, levels = c("CSF.1", "LIF.R", "CRP", "HbA1c-ebi-a-GCST90014006", "HbA1c-ieu-b-4842", "HDL-ieu-a-299", "HDL-ieu-b-109", "LDL-ebi-a-GCST90018961", "TG-ieu-a-302", "TG-ieu-b-111"))

plot_ind.effect <- indirect_effect_plot(data=all_traits, trait="Multimorbidity") # here plot product of coef 
plot_ind.effect + geom_point(shape = 1, size = 7.5,colour = "darkgrey") # final version of Figure 2 in the manuscript

# version 1 (with spelled out mediator names)
ggsave(paste0(PATH, "../mediation-product-coef-all-traits-superscript-", Sys.Date(), ".pdf"), height = 10, width = 7)
ggsave(paste0(PATH, "../mediation-product-coef-all-traits-superscript-", Sys.Date(), ".png"), height = 10, width = 7)

# version 2
# tidy mediator names
# all_traits$marker <- c("CSF.1", "LIF.R", "CRP", "HbA1c (GCST90014006)", "HbA1c (ieu-b-4842)", "HDL (ieu-a-299)", "HDL (ieu-b-109)", "LDL (GCST90018961)", "TG (ieu-a-302)", "TG (ieu-b-111)")
# # Change the levels of the marker variable to desired order in plot
# all_traits$marker <- factor(all_traits$marker, levels = c("CSF.1", "LIF.R", "CRP", "HbA1c (GCST90014006)", "HbA1c (ieu-b-4842)", "HDL (ieu-a-299)", "HDL (ieu-b-109)", "LDL (GCST90018961)", "TG (ieu-a-302)", "TG (ieu-b-111)" ))
# 
# plot_ind.effect.v2 <- indirect_effect_plot(data=all_traits, trait="Multimorbidity") # here plot product of coef 
# plot_ind.effect.v2                          

