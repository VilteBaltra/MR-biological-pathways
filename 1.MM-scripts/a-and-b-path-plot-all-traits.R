
#---------------------------------------------#
#                   SET UP
#---------------------------------------------#
rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

# specify path to dat
PATH_olink = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/olink/"
PATH_other = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/other-traits/"
PATH = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/"

# read in unadjusted a and b paths
data_a_b_paths_formatted_olink <- readRDS(paste0(PATH_olink, "data_a_b_paths_formatted_olink.rds") )
data_a_b_paths_formatted_other <- readRDS(paste0(PATH_other, "data_a_b_paths_formatted_other_traits.rds")) 

# subset olink to significant ones
data_a_b_paths_formatted_olink_sub <- data_a_b_paths_formatted_olink %>% filter(a.pval < 0.05 | b.pval < 0.05)

# combine output into one
data_a_b_paths_formatted <- rbind(data_a_b_paths_formatted_olink_sub, data_a_b_paths_formatted_other)

# specify names for each category
inflammatory_names <- c("CRP", data_a_b_paths_formatted_olink_sub$marker)
glycaemic_names <- c("Fasting glucose || id:ebi-a-GCST90002232","Fasting glucose || id:ieu-b-113","Fasting insulin || id:ebi-a-GCST90002238", "Fasting insulin || id:ieu-b-116","HbA1c || id:ebi-a-GCST90014006", "HbA1c || id:ieu-b-4842")
cortisol_name <- "Plasma cortisol || id:ieu-a-1012"
lipids_names <- as.character(unique((data_a_b_paths_formatted %>% filter(!marker %in% c(inflammatory_names, glycaemic_names, cortisol_name)))$marker))

# vector to specify the desired order
order_vector <- c(inflammatory_names, glycaemic_names, cortisol_name, lipids_names)

# arrange dataframe according to the order_vector
data_a_b_paths_formatted <- data_a_b_paths_formatted[match(order_vector, data_a_b_paths_formatted$marker), ]

# rename using shorten names
data_a_b_paths_formatted$marker <- c(inflammatory_names, "FG (GCST90002232)","FG (ieu-b-113)","FI (GCST90002238)", "FI (ieu-b-116)","HbA1c (GCST90014006)", "HbA1c (ieu-b-4842)",
                                     "Cortisol (ieu-a-1012)", "HDL (ieu-a-299)", "HDL (ieu-b-109)", "LDL (GCST90018961)", "LDL (ieu-a-300)", "TC (GCST90018974)", "TC (ieu-a-301)", "TG (ieu-a-302)", "TG (ieu-b-111)"    )

#---------------------------------#
#      PLOT a and b paths
#---------------------------------#
# stack data (for facet_wrap)
adata <- data_a_b_paths_formatted[, c("a", "a.se", "a.pval", "a.lCI", "a.hCI", "marker")]
adata$path <- "Maltreatment-Mediator"
colnames(adata) <- c("b", "b.se", "b.pval", "b.lCI", "b.hCI", "marker", "path") # make colnames same for rbind
bdata <- data_a_b_paths_formatted[, c("b", "b.se", "b.pval", "b.lCI", "b.hCI", "marker")]
bdata$path <- "Mediator-Multimorbidity"
abdata <- rbind(adata, bdata)

# # specify new names for each category
inflammatory_names # remains same
glycaemic_names_new <- c("FG (GCST90002232)","FG (ieu-b-113)","FI (GCST90002238)", "FI (ieu-b-116)","HbA1c (GCST90014006)", "HbA1c (ieu-b-4842)")
cortisol_name_new <- c("Cortisol (ieu-a-1012)")
lipids_names_new <- as.character(unique((abdata %>% filter(!marker %in% c(inflammatory_names, glycaemic_names_new, cortisol_name_new)))$marker))

# add category variable (used for assigning colours in plot)
abdata$category <- NA
abdata$category <- ifelse(abdata$marker %in% inflammatory_names, "Inflammatory",
                          ifelse(abdata$marker %in% glycaemic_names_new, "Glycaemic",
                                 ifelse(abdata$marker %in% cortisol_name_new, "Cortisol", "Lipids")))

# change the levels of the marker variable to desired order in plot
abdata$marker <- factor(abdata$marker, levels = c(inflammatory_names, cortisol_name_new, glycaemic_names_new, lipids_names_new) )

# save markers as factor and arrange accordingly
abdata$marker <- as.factor(abdata$marker) 
abdata <- abdata %>% arrange(marker) # makes sure order of pvalues and factor names are the same (important as the order of levels in the factor will determine the order in the plot)

# choose colours for forest plot
my_colors <- rep("#D9CDED", nrow(abdata)) # can use blue instead of purple or green "#a9def9"
# add colours for lipids
ind <- which(abdata$category == "Lipids")
my_colors[ind] <-"#EDE981" # "#fcf6bd" # lipid effects in yellow
# add colour for cortisol
ind <- which(abdata$category == "Cortisol")
my_colors[ind] <- "#d0f4de"  # "#EBCADF" # Cortisol effects in green
# add colours for glycaemic
ind <- which(abdata$category == "Glycaemic")
my_colors[ind] <- "#ff99c8" # "#EBCADF" # Glycaemic effects in pink


# reduce opacity of grid lines (default is 255)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 290)

# need to adjust so that a and b path significance is separated in plots
forest_ab_paths <- ggplot(
  data = abdata,
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
  coord_flip() +  # uncomment for vertical plot
  ggtitle(paste0(" ")) + # was Mediator-Multimorbidity
  theme(plot.title = element_text(size = 16, face = "bold")) + 
  scale_x_discrete(labels = c("C-reactive protein", inflammatory_names[2:length(inflammatory_names)], # spelling out CRP for clarity in figure
                              "Cortisol", 
                              expression("Fasting glucose"^c),expression("Fasting glucose"^b), expression("Fasting insulin"^c),expression("Fasting insulin"^b), expression("HbA1c"^c), expression("HbA1c"^b),
                              expression("HDL cholesterol"^a), expression("HDL cholesterol"^b), expression("LDL cholesterol"^c), expression("LDL cholesterol"^a), expression("Total cholesterol"^c), expression("Total cholesterol"^a), expression("Triglycerides"^a), expression("Triglycerides"^b)))

# search for "uncomment" keyword to know what to change for horizontal vs vertical plots

# main MS figure 1 (version 1)
temp <- forest_ab_paths +  facet_wrap(~ path, nrow = 1, scales="free_x", strip.position = "top") # uncomment for vertical plot
#facet_wrap(~ path, nrow = 2, scales="free_y")  # uncomment for horizontal plot

temp + geom_point(shape = 1,size = 8,colour = "darkgrey") # final version of Figure 1 in the manuscript
### save vertical ###
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-all-traits-superscript-", Sys.Date(), ".pdf"), height = 15, width = 15)
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-all-traits-superscript-", Sys.Date(), ".png"), height = 15, width = 15)

