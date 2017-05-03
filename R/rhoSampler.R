# Function to run the sampler through the rho parameter.
# It returns a list with three elements: MCMatrix (to update the matrix to the new values), currentValues (to update the vector of values to calculate the dendisties), nAccepted (to update the number of accepted values for each beta parameter),
# The arguments for the function are:
# data - A matrix with three columns. The first column contains the gene expression data for T1, the second contains the gene expression data for T2, and the third contains the eQTL values, V.
# MCMatrix - A matrix with the Markov chain along the columns and the parameters along the rows.
# currentValues - A vector that has the t - 1 values from MCMatrix. If a proposed value is accepted it is updated in the currentValues vector.
# proposalSD - A list containing values for the standard deviation for all possible parameters. In the case of the standard deviations the value represents the scale parameter.
# currentRho - An integer scalar used to select which rho parameter the sampler will run for.
# priorParameters - A list containing pairs of numbers for all possible parameters. The first number is the mean and the second number is the standard deviation. In the case of the priors for the standard deviations the two numbers are the shape and scale values.
# nAccepted - A matrix to hold the number of accepted values.
rhoSampler <- function (model,
                        data,
                        currentValues,
                        proposalSD,
                        currentRho,
                        priorParameters,
                        nAccepted) {

  # Create a vector to hold the prposed rho value
  newValues <- currentValues

  # Generate new value for current beta.
  rhoStar <- rtruncnorm(n = 1,
                        a = -1,
                        b = 1,
                        mean = currentValues[currentRho],
                        sd = proposalSD)

  newValues[currentRho] <- rhoStar

  # Calculate the numerator and denominator for the acceptance ratio
  numerator <- (log(dtruncnorm(x = newValues[currentRho],
                               a = -1,
                               b = 1,
                               mean = priorParameters[1],
                               sd = priorParameters[2]))

                + calcLogLikelihood(model = model,
                                    data = data,
                                    currentValues = newValues)

                + log(dtruncnorm(x = currentValues[currentRho],
                                 a = -1,
                                 b = 1,
                                 mean = newValues[currentRho],
                                 sd = proposalSD)))

  denominator <- (log(dtruncnorm(x = currentValues[currentRho],
                                 mean = priorParameters[1],
                                 sd = priorParameters[2]))

                  + calcLogLikelihood(model = model,
                                      data = data,
                                      currentValues = currentValues)

                  + log(dtruncnorm(x = newValues[currentRho],
                                   a = -1,
                                   b = 1,
                                   mean = currentValues[currentRho],
                                   sd = proposalSD)))

  acceptanceRatio <- numerator - denominator

  # Generate log uniform(0, 1) to compare to alpha which is min(acceptanceRatio, 0)
  logU <- log(runif(n = 1,
                    min = 0,
                    max = 1))

  alpha <- min(acceptanceRatio, 0)

  # Determine if the proposed beta should be accepted.
  if (!is.na(alpha) & logU < alpha) {

    currentValues <- newValues
    nAccepted <- nAccepted + 1

  }

  return (list(currentValues = currentValues,
               nAccepted = nAccepted))
}
