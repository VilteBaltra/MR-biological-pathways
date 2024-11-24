
#---------------------------------------------#
#                   SET UP
#---------------------------------------------#
rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

# specify path to dat
PATH_olink = " "
PATH_other = " "
PATH = " "
setwd(PATH)

# readPATH# read in unadjusted a and b paths
data_b_paths_formatted_olink <- openxlsx::read.xlsx("olink/unadjusted-b-path-md/olink-MD-uniMR-unadjusted-b-paths-2024-01-29.xlsx")
data_b_paths_formatted_other_old <- openxlsx::read.xlsx("other-traits/unadjusted-b-path-md/b-paths-unadjusted-Egger-Q-stat-uni-MR-2024-01-29.xlsx") 
data_b_paths_formatted_other <- openxlsx::read.xlsx("other-traits/unadjusted-b-path-md/lipids-large-with-and-withoutUKBB-DEPRESSION-uniMR-b-paths-2024-09-15.xlsx") 
data_b_crp <- openxlsx::read.xlsx("olink/unadjusted-b-path-md/CRP-MD-univariable-MR-unadjusted-b-2024-02-15.xlsx") 
data_b_metabolic <- openxlsx::read.xlsx("metabolites/combined-uniMR-b-paths-metabolites-DEPRESSION-2024-09-16.xlsx")

# keep ieu-b-4842 , ebi-a-GCST90014006, ieu-b-113,  ebi-a-GCST90002232, ebi-a-GCST90002238, ieu-b-116, ieu-a-1012
data_b_paths_formatted_other_old_sub <- data_b_paths_formatted_other_old %>% filter(id.exposure %in% c('ieu-b-4842', 'ebi-a-GCST90014006', 'ieu-b-113', 'ebi-a-GCST90002232', 'ebi-a-GCST90002238', 'ieu-b-116', 'ieu-a-1012'))
data_b_paths_formatted_other <- rbind(data_b_paths_formatted_other_old_sub, data_b_paths_formatted_other)  

# subset olink to significant ones
data_b_paths_formatted_olink$exposure <- sub("-1", "", data_b_paths_formatted_olink$exposure) # (remove '-1' from olink names)
data_b_paths_formatted_olink_sub <- data_b_paths_formatted_olink %>% filter(method == "Inverse variance weighted" & pval < 0.05)
data_b_paths_formatted_olink_sub <- data_b_paths_formatted_olink %>% filter(exposure %in% c(c('MIP.1.alpha', 'LAP.TGF.beta.1', 'HGF', 'FGF.21', 'CXCL6', 'CXCL5', 'CXCL1', 'CDCP1', 'CD244'), data_b_paths_formatted_olink_sub$exposure))

# subset metabolites to significant ones
data_b_metabolic_sub <- data_b_metabolic %>% filter(method == "Inverse variance weighted" & pval < 0.05)
data_b_metabolic_sub <- data_b_metabolic %>% filter(Exposure %in% c('Cit', data_b_metabolic_sub$Exposure))
data_b_metabolic_sub <- data_b_metabolic_sub[, c("Exposure", "Exposure.ID", "outcome", "method", "nsnp", "b", "se", "pval")] 
colnames(data_b_metabolic_sub) <- c("exposure", "id.exposure", "outcome", "method", "nsnp", "b", "se", "pval") 
data_b_metabolic_sub$id.outcome <- NA

# subset other traits to IVW
data_b_paths_formatted_other_sub <- data_b_paths_formatted_other %>% filter(method == "Inverse variance weighted")
data_b_crp <- data_b_crp[data_b_crp$method == "Inverse variance weighted",]

# combine output into one
data_b_paths_formatted <- rbind(data_b_paths_formatted_olink_sub, data_b_paths_formatted_other_sub, data_b_metabolic_sub, data_b_crp[, c('id.exposure', 'id.outcome', 'outcome', 'exposure', 'method', 'nsnp', 'b', 'se', 'pval')])

# specify names for each category
inflammatory_names <- c("CRP", unique(data_b_paths_formatted_olink_sub$exposure))
glycaemic_names <- c("Fasting glucose || id:ebi-a-GCST90002232","Fasting glucose || id:ieu-b-113","Fasting insulin || id:ebi-a-GCST90002238", "Fasting insulin || id:ieu-b-116", "Glycated haemoglobin HbA1c levels (UKB data field 30750) || id:ebi-a-GCST90014006", "HbA1c || id:ieu-b-4842")
cortisol_name <- "Cortisol || id:ieu-a-1012" 
lipids_names <- as.character(unique((data_b_paths_formatted %>% filter(!exposure %in% c(inflammatory_names, glycaemic_names, cortisol_name)))$exposure))

# vector to specify the desired order
order_vector <- c(inflammatory_names, glycaemic_names, cortisol_name, lipids_names)

# arrange dataframe according to the order_vector
data_b_paths_formatted <- data_b_paths_formatted[match(order_vector, data_b_paths_formatted$exposure), ]

#---------------------------------#
#      PLOT a and b paths
#---------------------------------#

# add category variable (used for assigning colours in plot)
data_b_paths_formatted$category <- NA
data_b_paths_formatted$marker <- data_b_paths_formatted$exposure
data_b_paths_formatted$category <- ifelse(data_b_paths_formatted$marker %in% inflammatory_names, "Inflammatory",
                          ifelse(data_b_paths_formatted$marker %in% glycaemic_names, "Glycaemic",
                                 ifelse(data_b_paths_formatted$marker %in% cortisol_name, "Cortisol", "Lipids")))

# change the levels of the marker variable to desired order in plot
data_b_paths_formatted$marker <- factor(data_b_paths_formatted$marker, levels = c(inflammatory_names, cortisol_name, glycaemic_names, lipids_names) )

# save markers as factor and arrange accordingly
data_b_paths_formatted$marker <- as.factor(data_b_paths_formatted$marker) 

# specify the order of levels of 'markers' as you wish to be plotted
data_b_paths_formatted$marker <- factor(data_b_paths_formatted$marker, 
                                        levels = c("CRP", inflammatory_names[2:length(inflammatory_names)], "Cortisol || id:ieu-a-1012", 
                                                   "bOHBut", "Cit", "DHA", "FAw3", "SM", "UnsatDeg",
                                                   "Fasting glucose || id:ebi-a-GCST90002232", "Fasting glucose || id:ieu-b-113", "Fasting insulin || id:ebi-a-GCST90002238",                                         
                                                   "Fasting insulin || id:ieu-b-116", "Glycated haemoglobin HbA1c levels (UKB data field 30750) || id:ebi-a-GCST90014006", 
                                                   "HbA1c || id:ieu-b-4842",  "nonHDL", "nonHDLwithUKBB","HDL", "HDLwithUKBB", "LDL", "LDLwithUKBB", "TC", "TCwithUKBB", "logTG", "logTGwithUKBB") )

#data_b_paths_formatted <- data_b_paths_formatted %>% arrange(marker) # makes sure order of pvalues and factor names are the same (important as the order of levels in the factor will determine the order in the plot)

# define data format function
my_data_format <- function(data){
  
  # add marker names
  #data$marker <- row.names(data)
  data$b.lCI <- data$b -1.96*data$se
  data$b.hCI <- data$b + 1.96*data$se
  
  # order alphabetically based on marker name
  #data <- data[order(data$marker),]
  return(data)
}

# format data for plotting
data_b_paths_formatted <- my_data_format(data=data_b_paths_formatted)

# choose colours for forest plot
my_colors <- rep("#D9CDED", nrow(data_b_paths_formatted)) # can use blue instead of purple or green "#a9def9"
# add colours for lipids
ind <- which(data_b_paths_formatted$category == "Lipids")
my_colors[ind] <-"#EDE981" # "#fcf6bd" # lipid effects in yellow
# add colour for cortisol
ind <- which(data_b_paths_formatted$category == "Cortisol")
my_colors[ind] <- "#d0f4de"  # "#EBCADF" # Cortisol effects in green
# add colours for glycaemic
ind <- which(data_b_paths_formatted$category == "Glycaemic")
my_colors[ind] <- "#EDE981" # was pink: "#ff99c8" # Glycaemic effects in pink


# reduce opacity of grid lines (default is 255)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 290)


# need to adjust so that a and b path significance is separated in plots
forest_b_paths <- ggplot(
  data = data_b_paths_formatted,
  aes(x = marker, y = b, ymin = b.lCI, ymax = b.hCI)
) +
  geom_errorbar(aes(ymin = b.lCI, ymax = b.hCI, col = marker),color = "darkgrey", width = 1, cex = 0.5) +
  geom_pointrange(aes(col = marker), size = 1.8,  color = my_colors) +
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed") +
  xlab(" ") + # was "b path"
  ylab("beta (95% CI)") +
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(color = NA, fill = NA),
    strip.text.y = element_text(face = "bold", size = 20),  # Increased font size here
    panel.grid.major.y = element_line(color = col_grid, size = 0.5), # was color = col_grid
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none",
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  # Rotate x-axis text vertically, uncomment for horizonal plot
    axis.text = element_text(face = "plain", size = 18),  # Increased font size for axis text
    axis.title = element_text(face = "plain", size = 20),  # Increased font size for axis title
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20)  # Increased font size for plot title
  ) + scale_color_manual(values = my_colors) + 
  coord_flip() + 
  ggtitle(paste0(" ")) + 
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(labels = c("C-reactive protein", inflammatory_names[2:length(inflammatory_names)], # spelling out CRP for clarity in figure
                              "Cortisol",
                              "3-hydroxybutyrate", "Citrate", "Docosahexaenoic acid", "Omega-3 fatty acids", "Sphingomyelins", "Degree of unsaturation",
                              expression("Fasting glucose"^a),expression("Fasting glucose"^b), expression("Fasting insulin"^a),expression("Fasting insulin"^b), expression("Glycated haemoglobin"^a), expression("Glycated haemoglobin"^b),
                              "non-HDL cholesterol (no UKBB)", "non-HDL cholesterol", "HDL cholesterol (no UKBB)", "HDL cholesterol", "LDL cholesterol (no UKBB)", "LDL cholesterol", 
                              "Total cholesterol (no UKBB)", "Total cholesterol", "Triglycerides (no UKBB)", "Triglycerides"))
# mediator-depression figure
forest_b_paths + geom_point(shape = 1,size = 8,colour = "darkgrey")

### save ###
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-all-traits-superscript-MD-", Sys.Date(), ".pdf"), height = 15, width = 10)
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-all-traits-superscript-MD-", Sys.Date(), ".png"), height = 15, width = 10)

