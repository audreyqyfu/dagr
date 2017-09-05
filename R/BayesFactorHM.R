# BayesFactorHM is a function that approximates Bayes Factor using the ratio of the harmonic mean of the posterior distribution of the likelihood between two models.
# It returns an approximation to Bayes Factor on either the log or normal scale
# the argument to the function are:
# firstModel - A string indicating which model will be in the numerator of Bayes Factor
# secondModel - A string indicating which model will be in the denominator of Bayes Factor
# firstMHOutput - The first element in the list of the output from the MHSampler function. This is the matrix containing the Markov chain for lambda.
# secondMHOutput - The first element in the list of the output from the MHSampler function. This is the matrix containing the Markov chain for lambda.
# data - A matrix with three columns. The first column contains the gene expression data for T1, the second contains the gene expression data for T2, and the third contains the genetic variant values, V.
# log - A logical argument indicating wether the results are reported on the log scale or the normal scale.
BayesFactorHM <- function (firstModel = 'model0',
                           secondModel = 'model1',
                           firstMHOutput,
                           secondMHOutput,
                           data,
                           log = FALSE) {

  # generate a value that is close to the log likelihood of the two models
  logA <- -(calcLogLikelihood(model = firstModel,
                              data = data,
                              currentValues = firstMHOutput) + 50)

  # Sum the likelihood for each column of the Markov chain matrix for the first model
  sumLikelihood1 <- 0
  for (e in 1:length(firstMHOutput[1, ])) {

    sumLikelihood1 <- sumLikelihood1 + (1 / exp(logA + calcLogLikelihood(model = firstModel,
                                                                         data = data,
                                                                         currentValues = firstMHOutput[, e])))

  }

  harmonicSum1 <- 1 / sumLikelihood1

  # Sum the likelihood for each column of the Markov chain matrix for the second model
  sumLikelihood2 <- 0
  for (v in 1:length(firstMHOutput[1, ])) {

    sumLikelihood2 <- sumLikelihood2 + (1 / exp(logA + calcLogLikelihood(model = secondModel,
                                                                         data = data,
                                                                         currentValues = secondMHOutput[, v])))

  }

  harmonicSum2 <- 1 / sumLikelihood2

  if (log == TRUE) {

    return (log(harmonicSum1) - log(harmonicSum2))

  } else {

    return (harmonicSum1 / harmonicSum2)

  }

}
