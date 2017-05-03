# TNoParent calculates the log likelihood when a gene has no parent
# It returns a scalar which is the log likelihood
# The arguments of the function are:
# geneExpression is a vector of gene expression data for the specified gene
# b0 is a scalar for the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
dTNoParent <- function(geneExpression, b0, sd) {
  return(dnorm(x = geneExpression,
               mean = b0,
               sd = sd,
               log = TRUE))
}




# TGivenV calculates the log likelihood when the only parent of a gene is V
# It returns a scalar which is the log likelihood
# The arguments of the function are:
# geneExpression is a vector of gene expression data for the specified gene
# V is a vector of eQTL values
# b0 is a scalar for the intercept of the linear model b0 + b1 * V which is the mean of the specified gene
# b1 is a scalar for the slope of the linear model b0 + b1 * V which is the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
dTGivenV <- function(geneExpression, V, b0, b1, sd) {
  return(dnorm(x = geneExpression,
               mean = b0 + b1 * V,
               sd = sd,
               log = TRUE))
}




# TGivenOtherT calculates the log likelihood when one gene is the parent of the other gene
# It returns a scalar which is the log likelihood
# The arguments of the function are:
# geneExpression1 is a vector of the gene expression for the parent gene
# geneExpression2 is a vector of the gene expression for the child gene
# b0 is a scalar for the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
# rho is the correlation of the two genes
# Mean is a scalar which is the emperical mean of the parent gene
# Sd is a scalar which is the emperical standard deviation of the parent gene
dTGivenOtherT <- function(geneExpression1, geneExpression2, b0, sd, rho, Mean, Sd) {
  return(dnorm(x = geneExpression2,
               mean = b0 + rho * (sd / Sd) * (geneExpression1 - Mean),
               sd = sd * sqrt(1 - rho^2),
               log = TRUE))
}




# TGivenVOtherT calculates the log likelihood when a gene has V and the other gene as parents
# It returns a scalar which is the log likelihood
# The arguments to the function are:
# geneExpression1 is a vector of the gene expression for the parent gene
# geneExpression2 is a vector of the gene expression for the child gene
# V is a vector of eQTL values
# b0 is a scalar for the intercept of the linear model b0 + b1 * V which is the mean of the specified gene
# b1 is a scalar for the slope of the linear model b0 + b1 * V which is the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
# rho is the correlation of the two genes
# Mean is a scalar which is the emperical mean of the parent gene
# Sd is a scalar which is the emperical standard deviation of the parent gene
dTGivenVOtherT <- function(geneExpression1, geneExpression2, V, b0, b1, sd, rho, Mean, Sd) {
  return(dnorm(x = geneExpression2,
               mean = b0 + b1 * V + rho * (sd / Sd) * (geneExpression1 - Mean),
               sd = sd * sqrt(1 - rho^2),
               log = TRUE))
}
