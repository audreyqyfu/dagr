# rTNoParent simulates N observations when a gene has no parent
# It returns a vector of length N
# The arguments of the function are:
# N is a scalar for the number of observations to simulate
# b0 is a scalar for the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
rTNoParent <- function (N,
                        b0,
                        sd) {

  return (rnorm(n = N,
                mean = b0,
                sd = sd))

}


# rTGivenV simulates N observations when the only parent of a gene is V
# It returns a vector of length N
# The arguments of the function are:
# N is a scalar for the number of observations to simulate
# V is a vector of eQTL values
# b0 is a scalar for the intercept of the linear model b0 + b1 * V which is the mean of the specified gene
# b1 is a scalar for the slope of the linear model b0 + b1 * V which is the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
rTGivenV <- function (N,
                      V,
                      b0,
                      b1,
                      sd) {

  return (rnorm(n = N,
                mean = b0 + b1 * V,
                sd = sd))

}


# rTGivenOtherT simulates N observations when one gene is the parent of the other gene
# It returns a vector of length N
# The arguments of the function are:
# N is a scalar for the number of observations to simulate
# geneExpression is a vector of the gene expression for the parent gene
# b0 is a scalar for the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
# rho is the correlation of the two genes
# Mean is a scalar which is the emperical mean of the parent gene
# Sd is a scalar which is the emperical standard deviation of the parent gene
rTGivenOtherT <- function (N,
                           geneExpression,
                           b0,
                           sd,
                           rho,
                           Mean,
                           Sd) {

  return (rnorm(n = N,
                mean = b0 + rho * (sd / Sd) * (geneExpression - Mean),
                sd = sd * sqrt(1 - rho^2)))

}


# rTGivenVOtherT simulates N observations when a gene has V and the other gene as parents
# It returns a vector of length N
# The arguments to the function are:
# N is a scalar for the number of observations to simulate
# V is a vector of eQTL values
# geneExpression is a vector of the gene expression for the parent gene
# b0 is a scalar for the intercept of the linear model b0 + b1 * V which is the mean of the specified gene
# b1 is a scalar for the slope of the linear model b0 + b1 * V which is the mean of the specified gene
# sd is a scalar for the standard deviation of the specified gene
# rho is the correlation of the two genes
# Mean is a scalar which is the emperical mean of the parent gene
# Sd is a scalar which is the emperical standard deviation of the parent gene
rTGivenVOtherT <- function (N, V, geneExpression, b0, b1, sd, rho, Mean, Sd) {

  return(rnorm(n = N,
               mean = b0 + b1 * V + rho * (sd / Sd) * (geneExpression - Mean),
               sd = sd * sqrt(1 - rho^2)))

}


# rTGivenVOtherTBivariate simulates N observations when a gene has V and the other gene as parents.
# It returns an N by 2 matrix.
# The arguments for the function are:
# b0.1 is a scalar for the intercept of the linear model b0.1 + b1.1 * V which is the mean of gene 1
# b0.2 is a scalar for the intercept of the linear model b0.2 + b1.2 * V which is the mean of gene 2
# b1.1 is a scalar for the slope of the linear model b0.1 + b1.1 * V which is the mean of gene 1
# b1.2 is a scalar for the slope of the linear model b0.2 + b1.2 * V which is the mean of gene 2
# sd.1 is a scalar for the standard deviation of gene 1
# sd.2 is a scalar for the standard deviation of gene 2
# rho is the correlation of the two genes
rTGivenVOtherTBivariate <- function (N,
                                     V,
                                     b0.1,
                                     b1.1,
                                     b0.2,
                                     b1.2,
                                     sd.1,
                                     sd.2,
                                     rho) {

  bivariateData <- matrix(nrow = N, ncol = 2)

  Sigma <- matrix(c(sd.1^2,
                    rho * sd.1 * sd.2,
                    rho * sd.1 * sd.2,
                    sd.2^2),
                  nrow = 2)

  for (e in 1:N) {

    bivariateData[e, ] <- rmvnorm(n = 1,
                                  mean = c(b0.1 + b1.1 * V[e],
                                           b0.2 + b1.2 * V[e]),
                                  sigma = Sigma)

  }

  return (bivariateData)

}
