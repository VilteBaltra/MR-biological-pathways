#---------------------------------------------#
#      CALCULATE INDIRECT EFFECT
#---------------------------------------------#
# Calculating indirect effect of all traits in olink-output-mediation folder 
library(tidyverse)
# define path to functions
FUNCTIONS = "~/Documents/Projects/MR-mediation-CM-MM/scripts/TIDY-WORKFLOW/source/" 
# path MM
PATH_2024 = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/other-traits/adjusted-b-path/"

# path to data 
BETA1 = paste0(PATH_2024, "selected-a-and-b-paths-adjusted-2024-02-10.xlsx")

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
openxlsx::write.xlsx(mediation, file = paste0("mediation-results-CM-other-traits-MM-adjusted-", Sys.Date(), ".xlsx"), rowNames = T, overwrite = T)
