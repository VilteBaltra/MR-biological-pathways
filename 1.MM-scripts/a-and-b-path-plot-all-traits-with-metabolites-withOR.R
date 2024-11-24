
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
PATH_metabolites = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/metabolites/"
PATH_other = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/other-traits/"
PATH = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/"

# read in unadjusted a and b paths
data_a_b_paths_formatted_olink <- readRDS(paste0(PATH_olink, "data_a_b_paths_formatted_olink.rds") )
data_a_b_paths_formatted_metabolites <- readRDS(paste0(PATH_metabolites, "data_a_b_paths_formatted_metabolites.rds") )
data_a_b_paths_formatted_other <- readRDS(paste0(PATH_other, "data_a_b_paths_formatted_other_traits_v2.rds")) 

# subset olink to significant ones
data_a_b_paths_formatted_olink_sub <- data_a_b_paths_formatted_olink %>% filter(a.pval < 0.05 | b.pval < 0.05)

# subset metabolites to significant ones
data_a_b_paths_formatted_metabolites_sub <- data_a_b_paths_formatted_metabolites %>% filter(a.pval < 0.05 | b.pval < 0.05)
#format to match other dfs
data_a_b_paths_formatted_metabolites_sub$outcome <- rownames(data_a_b_paths_formatted_metabolites_sub)
data_a_b_paths_formatted_metabolites_sub$marker <- rownames(data_a_b_paths_formatted_metabolites_sub)
data_a_b_paths_formatted_metabolites_sub$a.lCI <- data_a_b_paths_formatted_metabolites_sub$a - 1.96 * data_a_b_paths_formatted_metabolites_sub$a.se
data_a_b_paths_formatted_metabolites_sub$a.hCI <- data_a_b_paths_formatted_metabolites_sub$a + 1.96 * data_a_b_paths_formatted_metabolites_sub$a.se
data_a_b_paths_formatted_metabolites_sub$b.lCI <- data_a_b_paths_formatted_metabolites_sub$b - 1.96 * data_a_b_paths_formatted_metabolites_sub$b.se
data_a_b_paths_formatted_metabolites_sub$b.hCI <- data_a_b_paths_formatted_metabolites_sub$b + 1.96 * data_a_b_paths_formatted_metabolites_sub$b.se

# combine output into one
data_a_b_paths_formatted <- rbind(data_a_b_paths_formatted_olink_sub, data_a_b_paths_formatted_metabolites_sub[, colnames(data_a_b_paths_formatted_olink_sub)], data_a_b_paths_formatted_other)

# adding OR and CIs for b path only
# data_a_b_paths_formatted$a.OR <- exp(data_a_b_paths_formatted$a)
# data_a_b_paths_formatted$a.OR.lCI <- exp(data_a_b_paths_formatted$a.lCI)
# data_a_b_paths_formatted$a.OR.hCI <- exp(data_a_b_paths_formatted$a.hCI)

data_a_b_paths_formatted$b.OR <- exp(data_a_b_paths_formatted$b)
data_a_b_paths_formatted$b.OR.lCI <- exp(data_a_b_paths_formatted$b.lCI)
data_a_b_paths_formatted$b.OR.hCI <- exp(data_a_b_paths_formatted$b.hCI)

# specify names for each category
inflammatory_names <- c("CRP", data_a_b_paths_formatted_olink_sub$marker)
#glycaemic_names <- c("Fasting glucose || id:ebi-a-GCST90002232","Fasting glucose || id:ieu-b-113","Fasting insulin || id:ebi-a-GCST90002238", "Fasting insulin || id:ieu-b-116","HbA1c || id:ebi-a-GCST90014006", "HbA1c || id:ieu-b-4842")
cortisol_name <- "Plasma cortisol || id:ieu-a-1012"
metabolite_names <- as.character(unique((data_a_b_paths_formatted %>% filter(!marker %in% c(inflammatory_names, cortisol_name)))$marker))
metabolite_names <- c(metabolite_names[1:15], metabolite_names[22:31], metabolite_names[16:21])

# vector to specify the desired order
order_vector <- c(inflammatory_names, cortisol_name, metabolite_names)

# arrange dataframe according to the order_vector
data_a_b_paths_formatted <- data_a_b_paths_formatted[match(order_vector, data_a_b_paths_formatted$marker), ]

#---------------------------------#
#      PLOT a and b paths
#---------------------------------#
# stack data (for facet_wrap)
adata <- data_a_b_paths_formatted[, c("a", "a.se", "a.pval", "a.lCI", "a.hCI", "marker")]
adata$path <- "Maltreatment-Mediator"
colnames(adata) <- c("b", "b.se", "b.pval", "b.lCI", "b.hCI", "marker", "path") # make colnames same for rbind
bdata <- data_a_b_paths_formatted[, c("b.OR", "b.se", "b.pval", "b.OR.lCI", "b.OR.hCI", "marker")]
bdata$path <- "Mediator-Multimorbidity"
colnames(bdata) <- c("b", "b.se", "b.pval", "b.lCI", "b.hCI", "marker", "path") # make colnames same for rbind
abdata <- rbind(adata, bdata)

# # specify new names for each category
inflammatory_names # remains same

# add category variable (used for assigning colours in plot)
abdata$category <- NA
# abdata$category <- ifelse(abdata$marker %in% inflammatory_names, "Inflammatory",
#                           ifelse(abdata$marker %in% glycaemic_names_new, "Glycaemic",
#                                  ifelse(abdata$marker %in% cortisol_name_new, "Cortisol", "Lipids")))
abdata$category <- ifelse(abdata$marker %in% inflammatory_names, "Inflammatory",
                          ifelse(abdata$marker %in% cortisol_name, "Cortisol", "Metabolites"))

# change the levels of the marker variable to desired order in plot
abdata$marker <- factor(abdata$marker, levels = c(inflammatory_names, cortisol_name, metabolite_names) )

# save markers as factor and arrange accordingly
abdata$marker <- as.factor(abdata$marker) 
abdata <- abdata %>% arrange(marker) # makes sure order of pvalues and factor names are the same (important as the order of levels in the factor will determine the order in the plot)

# choose colours for forest plot
my_colors <- rep("#D9CDED", nrow(abdata)) # can use blue instead of purple or green "#a9def9"
# add colours for metabolites
ind <- which(abdata$category == "Metabolites")
my_colors[ind] <-"#EDE981" # "#fcf6bd" # lipid effects in yellow
# add colour for cortisol
ind <- which(abdata$category == "Cortisol")
my_colors[ind] <- "#d0f4de"  # "#EBCADF" # Cortisol effects in green

# reduce opacity of grid lines (default is 255)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 290)

# need to adjust so that a and b path significance is separated in plots
forest_ab_paths <- ggplot(
  data = abdata,
  aes(x = marker, y = b, ymin = b.lCI, ymax = b.hCI)
) +
  geom_errorbar(aes(ymin = b.lCI, ymax = b.hCI, col = marker),color = "darkgrey", width = 1, cex = 0.5) +
  geom_pointrange(aes(col = marker), size = 1.8,  color = my_colors) +
  geom_hline(data = abdata %>% filter(path == "Maltreatment-Mediator"), # new
             aes(yintercept = 0), color = "darkblue", linetype = "dashed") +
  geom_hline(data = abdata %>% filter(path == "Mediator-Multimorbidity"),
             aes(yintercept = 1), color = "darkblue", linetype = "dashed") + # new till here
  
  #geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed") +
  xlab(" ") + # was "b path"
  ylab(" ") +
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
  theme(plot.title = element_text(size = 16, face = "bold")) + # ***if any changes to input df below names (in scale_x_discrete) need to be updated accordingly!***
  scale_x_discrete(labels = c("C-reactive protein", inflammatory_names[2:length(inflammatory_names)], # spelling out CRP for clarity in figure
                              "Cortisol",
                              "Apolipoprotein A-I", "Apolipoprotein B",  "Citrate",  "Esterified cholesterol", "Omega-6 fatty acids",  "Free cholesterol",
                              "HDL particle diameter", "HDL2 cholesterol", "Linoleic acid", "LDL triglycerides", "Polyunsaturated fatty acids",
                              "Remnant cholesterol", "Sphingomyelins", "Total fatty acids", "VLDL cholesterol", # keeping nightingale metabolites the same
                              "HDL cholesterol (no UKBB)", "HDL cholesterol", "LDL cholesterol (no UKBB)", "LDL cholesterol",
                              "Tiglycerides (no UKBB)", "Triglycerides", "non-HDL cholesterol (no UKBB)", "non-HDL cholesterol", "Total cholesterol (no UKBB)", "Total cholesterol",
                              expression("Fasting glucose"^c),expression("Fasting glucose"^a), expression("Fasting insulin"^c),expression("Fasting insulin"^a), expression("Glycated haemoglobin"^a), expression("Glycated haemoglobin"^b)))

# main MS figure 1 (version 1)
# temp <- forest_ab_paths +  facet_wrap(~ path, nrow = 1, scales="free_x", strip.position = "bottom", labeller = as_labeller(c("Maltreatment-Mediator" = "beta (95% CI)", "Mediator-Multimorbidity" = "OR (95% CI)"))) # uncomment for vertical plot
# #facet_wrap(~ path, nrow = 2, scales="free_y")  # uncomment for horizontal plot
# 
# # Adjust the positioning of the facet labels to be under the numbers
temp <- forest_ab_paths + facet_wrap( ~ path, nrow = 1, scales = "free_x", strip.position = "bottom", 
    labeller = as_labeller(c(
      "Maltreatment-Mediator" = "beta (95% CI)", 
      "Mediator-Multimorbidity" = "OR (95% CI)"
    ))
  ) + 
  theme(
    strip.text = element_text(size = 20, face = "plain"),  # face = "bold" (for bold beta and or)
    strip.background = element_blank(),  # Remove background of the strip
    strip.placement = "outside",  # Move the strip outside
    panel.spacing = unit(1, "lines")  # Adjust spacing between panels if necessary
  )

# Generate the plot with the updated strip position
temp + geom_point(shape = 1, size = 8, colour = "darkgrey") +
  ggtitle("Maltreatment-Mediator                                             Mediator-Multimorbidity")

#temp + geom_point(shape = 1,size = 8,colour = "darkgrey") # final version of Figure 1 in the manuscript
### save vertical ###
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-all-traits-superscript-with-metabolites-OR-", Sys.Date(), ".pdf"), height = 19, width = 15) # with hbA1c spelled out
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-all-traits-superscript-with-metabolites-OR-", Sys.Date(), ".png"), height = 19, width = 15) # with hbA1c spelled out

