#---------------------------------------------#
#                   SET UP
#---------------------------------------------#
rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

# specify path to data
PATH = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/other-traits/"

# read in unadjusted a and b paths
a_df <- openxlsx::read.xlsx(paste0(PATH, "unadjusted-a-and-b-paths/a-and-b-paths-unadjusted-uni-MR-2024-01-31.xlsx"), sheet = 1, rowNames = F)
b_df_unadjusted <- openxlsx::read.xlsx(paste0(PATH, "unadjusted-a-and-b-paths/a-and-b-paths-unadjusted-uni-MR-2024-01-31.xlsx"), sheet = 2, rowNames = F)

# read in mediation results
mediation_results <- openxlsx::read.xlsx(paste0(PATH, "mediation-other-traits/mediation-results-CM-other-traits-MM-adjusted-2024-02-10.xlsx"), rowNames = T)
mediation_difference_method <- openxlsx::read.xlsx(paste0(PATH, "mediation-other-traits/mediation-results-difference-method-2024-01-24.xlsx"), rowNames = T)

#---------------------------------------------#
#                     PLOT
#---------------------------------------------#
# This script plots the MR mediation results (a path, b path and indirect effect)

# reduce opacity of grid lines (default is 255)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

# select only IVW methd for plotting
a2 <- a_df %>% filter(method == "Inverse variance weighted") 
# cortisol in b path has wald ratio, so rename to IVW to be included in plots
#b_df_unadjusted[71, "method"] <- "Inverse variance weighted" # ensure this is still 71 if rerunning
b2 <- b_df_unadjusted %>% filter(method == "Inverse variance weighted") 

# copy and rename key columns
a2[, c('a', 'a.se', 'a.pval')] <- a2[, c('b', 'se', 'pval')]
b2[, c('b', 'b.se', 'b.pval')] <- b2[, c('b', 'se', 'pval')]
# NOTE: double check 8th row is still triglycerides (also applied to below rows where trait names are shortened)
#b2[8, "exposure"] <- "Triglycerides || id:ieu-b-111" # change 'triglycerides' into capital T to match with a2 name

# combine output into same dataframe #  Triglycerides || id:ieu-a-302. Triglycerides || id:ieu-b-111
df_tidy <- merge(a2[, c('a', 'a.se', 'a.pval', 'outcome')], b2[, c('b', 'b.se', 'b.pval', 'exposure')], by.x='outcome', by.y='exposure')
df_tidy[6,"outcome"] <- "HbA1c || id:ebi-a-GCST90014006" # shortened name from "Glycated haemoglobin HbA1c levels (UKB data field 30750) || id:ebi-a-GCST90014006"
df_tidy[14, 'outcome'] <- "Total cholesterol || id:ebi-a-GCST90018974" # shortened name from 'Total cholesterol levels || id:ebi-a-GCST90018974'
rownames(df_tidy) <- df_tidy$outcome # add mediator name as row name

# define data format function
my_data_format <- function(data){

  # add marker names
  data$marker <- row.names(data)
  # add CIs for a path
  data$a.lCI <- data$a -1.96*data$a.se
  data$a.hCI <- data$a + 1.96*data$a.se
  # add CIs for b path
  data$b.lCI <- data$b -1.96*data$b.se
  data$b.hCI <- data$b + 1.96*data$b.se

  # order alphabetically based on marker name
  data <- data[order(data$marker),]
  return(data)
}

# format data
data_a_b_paths_formatted <- my_data_format(data=df_tidy)
#saveRDS(data_a_b_paths_formatted, file = paste0(PATH, "data_a_b_paths_formatted_other_traits.rds")) 

#---------------------------------#
#    PLOT a path (same for all)
#---------------------------------#

# identify significant effects for a path
ind <- which(data_a_b_paths_formatted$a.pval < 0.05)
# Choose colours for forest2 plot
my_colors <- rep("grey", nrow(data_a_b_paths_formatted))
my_colors[ind] <- "black" # significant effects in black

forest_plot_a_path <- ggplot(
  data = data_a_b_paths_formatted,
  aes(x = marker, y = a, ymin = a.lCI, ymax = a.hCI)
) +
  geom_pointrange(aes(col = marker)) +
  geom_hline(yintercept = 0, color = "red") +
  xlab("a path") +
  ylab("beta (95% CI)") +
  geom_errorbar(aes(ymin = a.lCI, ymax = a.hCI, col = marker), width = 0, cex = 0.5) +
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(color = NA, fill = NA),
    strip.text.y = element_text(face = "bold", size = 16),  # Increased font size here
    panel.grid.major.y = element_line(color = col_grid, size = 0.5),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none",
    axis.text = element_text(face = "plain", size = 14),  # Increased font size for axis text
    axis.title = element_text(face = "plain", size = 14),  # Increased font size for axis title
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18)  # Increased font size for plot title
  ) + scale_color_manual(values = my_colors) +
  coord_flip() + ggtitle("CM-Mediator") +
  theme(plot.title = element_text(size = 16, face = "bold"))  # Increased font size for plot title


#---------------------------------#
#           PLOT b path
#---------------------------------#

my_forest_plot_b_path <- function(data, trait) {
  # identify significant effects for b path
  ind <- which(data$b.pval < 0.05)
  # Choose colors for forest2 plot
  my_colors <- rep("grey", nrow(data))
  my_colors[ind] <- "black" # significant effects in black
  
  forest2 <- ggplot(
    data = data,
    aes(x = marker, y = b, ymin = b.lCI, ymax = b.hCI)
  ) +
    geom_pointrange(aes(col = marker)) +
    geom_hline(yintercept = 0, color = "red") +
    xlab("b path") +
    ylab("beta (95% CI)") +
    geom_errorbar(aes(ymin = b.lCI, ymax = b.hCI, col = marker), width = 0, cex = 0.5) +
    theme_classic() +
    theme(
      panel.background = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      strip.text.y = element_text(face = "bold", size = 16),  # Increased font size here
      panel.grid.major.y = element_line(color = col_grid, size = 0.5),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.position = "none",
      axis.text = element_text(face = "plain", size = 14),  # Increased font size for axis text
      axis.title = element_text(face = "plain", size = 14),  # Increased font size for axis title
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)  # Increased font size for plot title
    ) + scale_color_manual(values = my_colors) +
    coord_flip() + ggtitle(paste0("Mediator-", trait)) +
    theme(plot.title = element_text(size = 16, face = "bold"))  # Increased font size for plot title
}

# plot b path
forest_plot_b_path <- my_forest_plot_b_path(data_a_b_paths_formatted, "MM")

# combine a and b path plots into one
combined_plots <- forest_plot_a_path + forest_plot_b_path + plot_layout(ncol = 2)
print(combined_plots)

### save ###
ggsave(paste0(PATH, "../plots/a-and-b-paths-unadjusted-other-traits-", Sys.Date(), ".pdf"), height=10, width=15)

#---------------------------------#
#      PLOT indirect effect
#---------------------------------#
# # plot indirect effects from product of coef
# add same names as for a and b paths
mediation_results$marker <- row.names(mediation_results)
# mediation_results[2,"marker"] <- "Triglycerides || id:ieu-a-302" 
# mediation_results[3,"marker"] <- "HDL cholesterol || id:ieu-b-109"
# mediation_results[4,"marker"] <- "Triglycerides || id:ieu-b-111"
# mediation_results[5,"marker"] <- "HbA1c || id:ieu-b-4842"
# mediation_results[6,"marker"] <- "LDL cholesterol || id:ebi-a-GCST90018961"
# mediation_results[7,"marker"] <- "HbA1c || id:ebi-a-GCST90014006"
# row.names(mediation_results) <- mediation_results$marker

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

# define indirect effect plotting function
indirect_effect_plot <- function(data, trait, combined=FALSE) {
  # Identify significant effects for indirect effect
  ind <- which(data$ab.pval < 0.05)
  
  if (combined == TRUE) {
    
    # Choose colors for forest2 plot
    my_colors <- rep("black", nrow(data))
    # Create a mapping between unique names and numeric values
    unique_markers <- unique(data$marker)
    unique_numeric_values <- seq_along(unique_markers)
    mapping <- setNames(unique_numeric_values, unique_markers)
    # Create a new column 'marker_numeric' using the mapping
    data$marker_numeric <- mapping[data$marker]
  } else {
    # Choose colors for forest2 plot
    my_colors <- rep("grey", nrow(data))
    my_colors[ind] <- "black" # Significant effects in black
    data$marker_numeric <- data$marker
  }
  
  forest3 <- ggplot(
    data = subset(data, difference.method == "NO"),
    aes(x = factor(marker_numeric), y = ab, ymin = ab.lCI, ymax = ab.hCI)
  ) +
    geom_pointrange(
      # data = subset(data, difference.method == "NO"),
      aes(col = marker)) +
    geom_hline(yintercept = 0, color = "red") +
    xlab("Indirect effect") +
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
    ) + scale_color_manual(values = my_colors) +
    coord_flip() + 
    ggtitle(paste0("CM-", trait)) +
    theme(plot.title = element_text(size = 16, face = "bold"))
  
  if (combined == TRUE) {
    forest3 <- forest3 +
      geom_pointrange(
        data = subset(data, difference.method == "YES"),
        aes(x = marker_numeric-0.2, y = ab, ymin = ab.lCI, ymax = ab.hCI),
        col = "pink"
      ) + scale_x_discrete(labels = unique_markers) # or data$markers
  }
  return(forest3)
}


# to plot indirect effect (product of coef together with the difference method)
# first name markers equally to 'mediation_results_formatted' (make sure order of markers same if rerunning)
mediation_difference_method$marker <- c('CRP', 'HbA1c || id:ieu-b-4842', 'HDL cholesterol || id:ieu-b-109', 'Triglycerides || id:ieu-a-302', 'Triglycerides || id:ieu-b-111')
# first add same column names to 'mediation_difference_method' as in 'mediation_results_formatted' 
mediation_difference_method[, c('ab', 'ab.lCI', 'ab.hCI', 'ab.pval', 'marker')] <- mediation_difference_method[, c('indirect.effect', 'lCI', 'hCI', 'p_value', 'marker')]
# then add a column to distinguish between difference and product of coef methods
mediation_difference_method$difference.method <- "YES"
mediation_results_formatted$difference.method <- "NO"
# now combine the two dataframes 
combined_med_df <- rbind(mediation_results_formatted[, c('ab', 'ab.lCI', 'ab.hCI', 'ab.pval', 'marker', 'difference.method')], mediation_difference_method[, c('ab', 'ab.lCI', 'ab.hCI', 'ab.pval', 'marker', 'difference.method')]) 
print(combined_med_df)

## save ##
# obtain and save the plots (combined = TRUE / FALSE indicates whether product of coef should be plotted together with difference method)
plot1 <- indirect_effect_plot(data=combined_med_df, trait="Multimorbidity", combined = FALSE) # here only plot product of coef estimates
ggsave(paste0(PATH, "../plots/mediation-product-coef-other-traits-", Sys.Date(), ".pdf"), height=8, width=7.5)
#plot2 <- indirect_effect_plot(data=combined_med_df, trait="Multimorbidity", combined = TRUE) # here plot both product of coef and difference method estimates
#ggsave(paste0(PATH, "../plots/mediation-both-methods-other-traits-", Sys.Date(), ".pdf"), height=8, width=11)

