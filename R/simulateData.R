# simulateData generates data for gene 1, gene 2, and the genotype.
# It returns a dataframe with three columns. The first two columns are the gene expression data for the two genes and the third column is the genotype data.
# The arguments for the function are:
# N - A scalar for the number of observations to generate.
# p - A real number between 0 and 1 which represents the proportion of the population that has the reference allele
# seed - A positive integer for random generation
# model - A character string specifying which model to generate data for. It is of the form 'modelx' where x is a number between 0 and 4.
# parameters - A list of parameters used to generate the data (b0.1, b1.1, b0.2, b1.2, sd.1, sd.2, rho)
# For every model b0.1, b1.1, and sd.1 are the parameters to simulate data for T1
# Likewise b0.2, b1.2, and sd.2 are the parameters to simulate data for T2
simulateData <- function (N,
                          p,
                          seed = 338,
                          model = 'model0',
                          parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0)) {

  set.seed(seed)

  # V is a vector representing the genotype of an individual.
  # It follows a multinomial distribution with the following probabilities:
  # 0 - (1 - p)^2
  # 1 - 2 * p * (1 - p)
  # 2 - p^2
  # Where p is the probability of seeing the reference allele
  V <- c(sample(c(0, 1, 2),
                size = N,
                replace = TRUE,
                prob = c((1 - p)^2,
                         2*p*(1 - p),
                         p^2)))

  # The switch function is used to select which model to simulate data for based on the input of the argument, model, in the simulateData function.
  # It returns a dataframe with the expression level for gene 1 in the first column, the expression level for gene 2 in the second column, and the genotype in the third column.
  switch(model,

         # The mean for T1 given V is b0.1 + b1.1 * V
         # The standard deviation for T1 given V is sd.1
         # The mean for T2 is b0.2
         # The standard deviation for T2 is sd.2
         model0 = {

           T1 <- rTGivenV(N = N,
                          V = V,
                          b0 = parameters$b0.1,
                          b1 = parameters$b1.1,
                          sd = parameters$sd.1)

           T2 <- rTNoParent(N = N,
                            b0 = parameters$b0.2,
                            sd = parameters$sd.2)

           return (data.frame(T1.0 = T1,
                              T2.0 = T2,
                              V))

         },


         # The mean for T1 given V is b0.1 + b1.1 * V
         # The standard deviation for T1 given V is sd.1
         # The mean for T2 given T1 is b0.2 + rho * (sd.2 / sd(T1)) * (t1 - mean(T1))
         # The standard deviation for T2 given T1 is sd.2 * sqrt(1 - rho^2)
         model1 = {

           T1 <- rTGivenV(N = N,
                          V = V,
                          b0 = parameters$b0.1,
                          b1 = parameters$b1.1,
                          sd = parameters$sd.1)

           T2 <- rTGivenOtherT(N = N,
                               geneExpression = T1,
                               b0 = parameters$b0.2,
                               sd = parameters$sd.2,
                               rho = parameters$rho,
                               Mean = mean(T1),
                               Sd = sd(T1))

           return (data.frame(T1.1 = T1,
                              T2.1 = T2,
                              V))

         },


         # The mean for T1 given V and T2 is b0.1 + b1.1 * V + rho * (sd.1 / sd(T2)) * (t2 - mean(T2))
         # The standard deviation for T1 given V and T2 is sd.1 * sqrt(1 - rho^2)
         # The mean for T2 is b0.2
         # The standard deviation for T2 is sd.2
         model2 = {

           T2 <- rTNoParent(N = N,
                            b0 = parameters$b0.2,
                            sd = parameters$sd.2)

           T1 <- rTGivenVOtherT(N = N,
                                V = V,
                                geneExpression = T2,
                                b0 = parameters$b0.1,
                                b1 = parameters$b1.1,
                                sd = parameters$sd.1,
                                rho = parameters$rho,
                                Mean = mean(T2),
                                Sd = sd(T2))

           return (data.frame(T1.2 = T1,
                              T2.2 = T2,
                              V))

         },


         # The mean for T1 given V is b0.1 + b1.1 * V
         # The standard deviation for T1 given V is sd.1
         # The mean for T2 given V is b0.2 + b1.2 * V
         # The standard deviation for T2 given V is sd.2
         model3 = {

           T1 <- rTGivenV(N = N,
                          V = V,
                          b0 = parameters$b0.1,
                          b1 = parameters$b1.1,
                          sd = parameters$sd.1)

           T2 <- rTGivenV(N = N,
                          V = V,
                          b0 = parameters$b0.2,
                          b1 = parameters$b1.2,
                          sd = parameters$sd.2)

           return (data.frame(T1.3 = T1,
                              T2.3 = T2,
                              V))

         },


         # The mean for T1 given V is b0.1 + b1.1 * V
         # The mean for T2 given V is b0.2 + b1.2 * V
         # The covariance matrix is | sd.1^2, rho * sd.1 * sd.2 |
         #                          | rho * sd.1 * sd.2, sd.2^2 |
         model4 = {

           T1T2 <- rTGivenVOtherTBivariate(N = N,
                                           V = V,
                                           b0.1 = parameters$b0.1,
                                           b1.1 = parameters$b1.1,
                                           b0.2 = parameters$b0.2,
                                           b1.2 = parameters$b1.2,
                                           sd.1 = parameters$sd.1,
                                           sd.2 = parameters$sd.2,
                                           rho = parameters$rho)

           return (data.frame(T1.4 = T1T2[, 1],
                              T2.4 = T1T2[, 2],
                              V))

         },

         stop ("Model not included or missing"))
}
