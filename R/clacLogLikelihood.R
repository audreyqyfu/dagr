# calcLogLikelihood calculates the log likelihood for a specified model
# It returns a scalar which is the log likelihood of T1 and T2 added together
# The arguments for the function are:
# model - A character string specifiying which model to compute the log likelihood for# data is a matrix of three columns. The first two columns are gene expression data. The first column is data for gene 1, the second column is data for gene 2 and the third column has the eQTL values
# data - A matrix with three columns. The first column contains the gene expression data for T1, the second contains the gene expression data for T2, and the third contains the eQTL values, V.
# currentValues - A list of parameters (b0.1, b1.1, b0.2, b1.2, sd.1, sd.2, rho). the list for currentValues will change depending on which model we are running the sampler for.
calcLogLikelihood <- function (model,
                               data,
                               currentValues) {

  switch (model,

          # The mean for T1 given V is b0.1 + b1.1 * V
          # The standard deviation for T1 given V is sd.1
          # The mean for T2 is b0.2
          # The standard deviation for T2 is sd.2
          model0 = {

            T1 <- dTGivenV(geneExpression = data[, 1],
                           V = data[, 3],
                           b0 = currentValues[1],
                           b1 = currentValues[2],
                           sd = currentValues[4])

            T2 <- dTNoParent(geneExpression = data[, 2],
                             b0 = currentValues[3],
                             sd = currentValues[5])

            return (sum(T1 + T2))

          },

          # The mean for T1 given V is b0.1 + b1.1 * V
          # The standard deviation for T1 given V is sd.1
          # The mean for T2 given T1 is b0.2 + rho * (sd.2 / sd(T1)) * (t1 - mean(T1))
          # The standard deviation for T2 given T1 is sd.2 * sqrt(1 - rho^2)
          model1 = {

            T1 <- dTGivenV(geneExpression = data[, 1],
                           V = data[, 3],
                           b0 = currentValues[1],
                           b1 = currentValues[2],
                           sd = currentValues[4])

            Mean <- mean(data[, 1])
            Sd <- sd(data[, 1])

            T2 <- dTGivenOtherT(geneExpression1 = data[, 1],
                                geneExpression2 = data[, 2],
                                b0 = currentValues[3],
                                sd = currentValues[5],
                                rho = currentValues[6],
                                Mean = Mean,
                                Sd = Sd)

            return (sum(T1 + T2))

          },

          # The mean for T1 given V and T2 is b0.1 + b1.1 * V + rho * (sd.1 / sd(T2)) * (t2 - mean(T2))
          # The standard deviation for T1 given V and T2 is sd.1 * sqrt(1 - rho^2)
          # The mean for T2 is b0.2
          # The standard deviation for T2 is sd.2
          model2 = {

            Mean <- mean(data[, 2])
            Sd <- sd(data[, 2])

            T1 <- dTGivenVOtherT(geneExpression1 = data[, 2],
                                 geneExpression2 = data[, 1],
                                 V = data[, 3],
                                 b0 = currentValues[1],
                                 b1 = currentValues[2],
                                 sd = currentValues[4],
                                 rho = currentValues[6],
                                 Mean = Mean,
                                 Sd = Sd)

            T2 <- dTNoParent(geneExpression = data[, 2],
                             b0 = currentValues[3],
                             sd = currentValues[5])

            return (sum(T1 + T2))

          },

          # The mean for T1 given V is b0.1 + b1.1 * V
          # The standard deviation for T1 given V is sd.1
          # The mean for T2 given V is b0.2 + b1.2 * V
          # The standard deviation for T2 given V is sd.2
          model3 = {

            T1 <- dTGivenV(geneExpression = data[, 1],
                           V = data[, 3],
                           b0 = currentValues[1],
                           b1 = currentValues[2],
                           sd = currentValues[5])

            T2 <- dTGivenV(geneExpression = data[, 2],
                           V = data[, 3],
                           b0 = currentValues[3],
                           b1 = currentValues[4],
                           sd = currentValues[6])

            return (sum(T1 + T2))

          },

          # The mean for T1 given V is b0.1 + b1.1 * V
          # The standard deviation for T1 given V is sd.1
          # The mean for T2 given V and T1 is b0.2 + b1.2 * V + rho * (sd.2 / sd(T1)) * (t1 - mean(T1))
          # The standard deviation for T2 given V and T1 is sd.2 * sqrt(1 - rho^2)
          model4 = {

            T1T2 <- bivariateNormalDensity(data = data,
                                           currentValues = currentValues)

            return (sum(T1T2))

          }
  )
}
