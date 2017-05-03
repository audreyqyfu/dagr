# bivariateNormalDensity calculates the log likelihood for a specified model
# It returns a scalar which is the log likelihood of the bivariate distribution of T1 and T2 added together
# The arguments for the function are:
# data - A matrix with three columns. The first column contains the gene expression data for T1, the second contains the gene expression data for T2, and the third contains the eQTL values, V.
# currentValues - A list of parameters (b0.1, b1.1, b0.2, b1.2, sd.1, sd.2, rho)
bivariateNormalDensity <- function (data,
                                    currentValues) {

  x1 <- data[, 1]
  x2 <- data[, 2]
  mu1 <- currentValues[1] + currentValues[2] * data[, 3]
  mu2 <- currentValues[3] + currentValues[4] * data[, 3]
  sigma1 <- currentValues[5]
  sigma2 <- currentValues[6]
  rho <- currentValues[7]

  z <- (x1 - mu1)^2 / sigma1^2 - (2 * rho * (x1 - mu1) * (x2 - mu2)) / (sigma1 * sigma2) + (x2 - mu2)^2 / sigma2^2
  probDensityFunction <- -z / (2 * (1 - rho^2)) - log(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2))

  return (probDensityFunction)

}
