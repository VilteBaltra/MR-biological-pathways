# function to calculate the indirect effect, SE and p-value, using the product of coefficients method

product.coef <- function(a, b, a.se, b.se, a.pval, b.pval, total.effect){
  
  # LDL product of coef method calculation
  ab = a*b
  # ab.se = sqrt((SEa^2 * b^2) + (SEb^2 * a^2))
  ab.se = sqrt((a.se^2 * b^2) + (b.se^2 * a^2))
  z = ab / ab.se
  #ab.pval = exp(-0.717*z - 0.416*z^2)
  ab.pval = 2*pnorm(q=z, lower.tail=FALSE)
  # proportion mediated  
  prop = ab / total.effect * 100 
  # add CIs
  # To construct 95% confidence intervals for β, we take the coefficient and 
  # add/subtract 1.96 × SE of β. This is because β is assumed to follow a normal distribution 
  # and 95% of the values in the sampling distribution are contained within 1.96 standard errors of β.
  lCI = round((ab - 1.96 * ab.se), 3)
  hCI = round((ab + 1.96 * ab.se), 3)
  OR = round(exp(ab), 3)
  ORlCI = round(exp(ab - 1.96 * ab.se), 3)
  ORhCI = round(exp(ab + 1.96 * ab.se), 3)
  # combine all info
  out = round(cbind(a, a.se, a.pval, b, b.se, b.pval, ab, ab.se, lCI, hCI, ab.pval, prop, OR, ORlCI, ORhCI), 5)
  return(as.data.frame(out))
  
}
