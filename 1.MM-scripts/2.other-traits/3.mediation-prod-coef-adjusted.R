#---------------------------------------------#
#      CALCULATE INDIRECT EFFECT
#---------------------------------------------#
# Calculating indirect effect of all traits in olink-output-mediation folder 
library(tidyverse)
# define path to functions
FUNCTIONS = " " 
# path MM
PATH_2024 = " "

# path to data 
BETA1 = paste0(PATH_2024, "selected-a-and-b-paths-adjusted-2024-09-13.xlsx")

# read in a and b path data 
a_and_b <- openxlsx::read.xlsx(BETA1, sheet = 1, rowNames=T)

# source product of coefficients function
source(paste0(FUNCTIONS, "product.coef.R")) 

mediation <- product.coef(a = a_and_b$a, b = a_and_b$b, 
                          a.se=a_and_b$a.se, b.se=a_and_b$b.se, 
                          a.pval=a_and_b$a.pval, b.pval=a_and_b$b.pval,
                          total.effect = 0.15975) # total effect for CM-MM added
rownames(mediation) <- rownames(a_and_b)
print(mediation) # print all output
mediation[mediation$ab.pval < 0.05,] # print only significant mediators

# save output
setwd(PATH_2024)
openxlsx::write.xlsx(mediation, file = paste0("../mediation-other-traits/mediation-results-CM-other-traits-MM-adjusted-", Sys.Date(), ".xlsx"), rowNames = T, overwrite = T)
