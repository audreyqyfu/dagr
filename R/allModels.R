# allModels generates data for gene 1 and gene 2 under all five models. It also generates data for the genotype. The same vector for the genotype is used in generating the expression level data for all five models.
# It returns a dataframe with 11 columns. The first ten columns are the gene expression levels for the five models and the eleventh column is the genotype.
# The arguments for the function are:
# N - A scalar for the number of observations to generate.
# p - A real number between 0 and 1 which represents the proportion of the population that has the reference allele
# seed - A positive integer for random generation
# parameters - A list of parameters used to generate the data (b0.1, b1.1, b0.2, b1.2, sd.1, sd.2, rho)
allModels <- function (N, p, seed, parameters) {

  Models <- list()

  for (e in 1:5) {

    Models[[e]] <- simulateData(N = N,
                                p = p,
                                seed = seed,
                                paste('model', e - 1, sep = ''),
                                parameters = parameters)
  }

  return (cbind(Models[[1]][, 1:2], Models[[2]][, 1:2], Models[[3]][, 1:2], Models[[4]][, 1:2], Models[[5]]))
}
