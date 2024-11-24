#---------------------------------------------#
#                   SET UP
#---------------------------------------------#
rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

# PATH
PATH <- "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/metabolites/"

# read in unadjusted a and b paths
a_df <- openxlsx::read.xlsx(paste0(PATH, "unadjusted-a-and-b-path-metabolites/metabolites-a-paths-unadjusted-uni-MR-2024-08-09-withnames-nothirdsource.xlsx"), sheet = 1, rowNames = F)
b_df_unadjusted <- openxlsx::read.xlsx(paste0(PATH, "unadjusted-a-and-b-path-metabolites/metabolites-b-paths-unadjusted-uni-MR-2024-08-09-nothirdsource.xlsx"), sheet = 1, rowNames = F)

#---------------------------------------------#
#                     FORMAT
#---------------------------------------------#
# This script plots the MR mediation results (a path, b path and indirect effect)

# reduce opacity of grid lines (default is 255)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

# select only IVW methd for plotting
a2 <- a_df %>% filter(method == "Inverse variance weighted") 
b2 <- b_df_unadjusted %>% filter(method == "Inverse variance weighted") 

# copy and rename key columns
a2[, c('a', 'a.se', 'a.pval')] <- a2[, c('b', 'se', 'pval')]
b2[, c('b', 'b.se', 'b.pval')] <- b2[, c('b', 'se', 'pval')]

# combine output into same dataframe 
df_tidy <- merge(a2[, c('a', 'a.se', 'a.pval', 'outcome')], b2[, c('b', 'b.se', 'b.pval', 'exposure')], by.x='outcome', by.y='exposure')
#rownames(df_tidy) <- sub("-1", "", df_tidy$outcome) # add mediator name as row name (remove '-1' from marker names)

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
data_a_b_paths_formatted$marker <- data_a_b_paths_formatted$outcome
#saveRDS(data_a_b_paths_formatted, file = paste0(PATH, "data_a_b_paths_formatted_olink.rds")) 

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

# arrange according to marker
abdata <- abdata %>% arrange(marker) # makes sure order of pvalues and factor names are the same (important as the order of levels in the factor will determine the order in the plot)

# identify significant effects for a path
ind <- which(abdata$b.pval < 0.05)
# Choose colours for forest2 plot
my_colors <- rep("#D9CDED", nrow(data_a_b_paths_formatted))

# need to adjust so that a and b path significance is separated in plots
forest_ab_paths <- ggplot(
  data = abdata,
  aes(x = marker, y = b, ymin = b.lCI, ymax = b.hCI)
) +
  geom_pointrange(aes(col = marker), size = 1) +
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed") +
  xlab(" ") + # was "b path"
  ylab("beta (95% CI)") +
  geom_errorbar(aes(ymin = b.lCI, ymax = b.hCI, col = marker), width = 1, cex = 0.5) +
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(color = NA, fill = NA),
    strip.text.y = element_text(face = "bold", size = 16),  # Increased font size here
    panel.grid.major.y = element_line(color = col_grid, size = 0.5),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none",
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  # Rotate x-axis text vertically, uncomment for horizonal plot
    axis.text = element_text(face = "plain", size = 14),  # Increased font size for axis text
    axis.title = element_text(face = "plain", size = 14),  # Increased font size for axis title
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18)  # Increased font size for plot title
  ) + scale_color_manual(values = my_colors) + 
  coord_flip() +  # uncomment for vertical plot
  ggtitle(paste0(" ")) + # was Mediator-Multimorbidity
  theme(plot.title = element_text(size = 16, face = "bold")) + 
  #facet_wrap(~ path, nrow = 2, scales="free_y")  # uncomment for horizontal plot
  facet_wrap(~ path, nrow = 1, scales="free_x", strip.position = "top") # uncomment for vertical plot

forest_ab_paths # search for "uncomment" keyword to know what to change for horizontal vs vertical plots


# Apply specific color scale to path A facet
abdata_path1 <- subset(abdata, path == "Maltreatment-Mediator")
ind <- which(abdata_path1$b.pval < 0.05)
my_colors3 <- rep("#D9CDED", nrow(abdata_path1))
my_colors3[ind] <- "#845EC2" # significant effects in darker purple
forest_ab_paths2 <- forest_ab_paths +
  geom_pointrange(
    data = abdata_path1,
    aes(col = marker), size = 1,
    color = my_colors3
  )

# Apply specific color scale to path B facet
abdata_path2 <- subset(abdata, path == "Mediator-Multimorbidity")
ind <- which(abdata_path2$b.pval < 0.05)
my_colors3 <- rep("#D9CDED", nrow(abdata_path2))
my_colors3[ind] <- "#845EC2" # significant effects in darker purple
forest_ab_paths3 <- forest_ab_paths2 +
  geom_pointrange(
    data = abdata_path2,
    aes(col = marker), size = 1,
    color = my_colors3
  )
# view plot
forest_ab_paths3

### save vertical ###
ggsave(paste0(PATH, "a-and-b-paths-unadjusted-metabolites-", Sys.Date(), ".pdf"), height = 18, width = 12)

