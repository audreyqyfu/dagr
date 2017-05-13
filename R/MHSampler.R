# MHSampler is the main function that will run the Metropolis Hastings sampler for a specified model.
# It returns a list. The first element of the list is a matrix with a Markov Chain for each parameter in the model. The second element is a matrix with the number of accepted values for each parameter.
# The arguments for the function are:
# nIterations - A scalar indicating the number of times the sampler runs.
# model - A character string specifying which model to run the sampler for.
# data - A matrix with three columns. The first column contains the gene expression data for T1, the second contains the gene expression data for T2, and the third contains the eQTL values, V.
# parameters - A list with all possible parameters. The values entered serve as the starting point for the sampler
# priorParameters - A list containing pairs of numbers for all possible parameters. The first number is the mean and the second number is the standard deviation. In the case of the priors for the standard deviations the two numbers are the shape and scale values.
# proposalSD - a list containing values for the standard deviation for all possible parameters. In the case of the standard deviations the value represents the scale parameter.
MHSampler <- function (nIterations = 100,
                       burnIn = 0.2,
                       thinned = 300,
                       model = 'model0',
                       data,
                       parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72),
                       priorParameters = list(b0.1 = c(0, 1), b1.1 = c(0, 1), b0.2 = c(0, 1), b1.2 = c(0, 1), sd.1 = c(1, 1), sd.2 = c(1, 1), rho = c(0, 0.5)),
                       proposalSD = list(b0.1 = 1, b1.1 = 1, b0.2 = 1, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.5)) {

  if (model == 'model1' | model == 'model2') {

    # Vector of places that are present in a specific model.
    parameterPlaces <- c(1:3, 5:7)

    # Create indicies to change the start and stop values of the for loops for the betaSampler, sigmaSamper, and rhoSampler functions.
    endBeta <- 3
    startSigma <- 4
    endSigma <- 5
    startRho <- 6

  } else if (model == 'model0') {

    # Vector of places that are present in a specific model.
    parameterPlaces <- c(1:3, 5:6)

    # Create indicies to change the start and stop values of the for loops for the betaSampler and sigmaSampler functions.
    endBeta <- 3
    startSigma <- 4
    endSigma <- 5

  } else if (model == 'model3') {

    # Vector of places that are present in a specific model.
    parameterPlaces <- c(1:6)

    # Create indicies to change the start and stop values of the for loops for the betaSampler, sigmaSamper, and rhoSampler functions.
    endBeta <- 4
    startSigma <- 5
    endSigma <- 6

  } else {

    # Vector of places that are present in a specific model.
    parameterPlaces <- c(1:7)

    # Create indicies to change the start and stop values of the for loops for the betaSampler, sigmaSamper, and rhoSampler functions.
    endBeta <- 4
    startSigma <- 5
    endSigma <- 6
    startRho <- 7

  }

  # Change the length of parameters, priorParameters, and proposalSD to match the number parameters used in model 1
  parameters = parameters[parameterPlaces]
  priorParameters =  priorParameters[parameterPlaces]
  proposalSD = proposalSD[parameterPlaces]

  # Matrix to hold the sample for each parameter
  MCMatrix <- matrix(nrow = length(parameters), ncol = nIterations)
  MCMatrix[, 1] <- unlist(parameters, use.names = FALSE)

  # Vector to hold the number of accepted values
  nAccepted <- c(rep(0, length(parameters)))

  for (t in 2:nIterations) {

    # Create a vector to hold the previously accepted values and to take the proposed values that are accepted
    currentValues <- MCMatrix[, t - 1]

    for (i in 1:endBeta) {
      proposedBeta <- betaSampler(model = model,
                                  data = data,
                                  currentValues = currentValues,
                                  proposalSD = proposalSD[[i]],
                                  currentBeta = i,
                                  priorParameters = priorParameters[[i]],
                                  nAccepted = nAccepted[i])

      # Update the currentValues vector and the nAccepted vector with the new quantities
      currentValues <- proposedBeta$currentValues
      nAccepted[i] <- proposedBeta$nAccepted
    }

    for (i in startSigma:endSigma) {
      proposedSigma <- sigmaSampler(model = model,
                                    data = data,
                                    currentValues = currentValues,
                                    proposalSD = proposalSD[[i]],
                                    currentSigma = i,
                                    priorParameters = priorParameters[[i]],
                                    nAccepted = nAccepted[i])

      # Update the currentValues vector and the nAccepted vector with the new quantities
      currentValues <- proposedSigma$currentValues
      nAccepted[i] <- proposedSigma$nAccepted
    }

    if (model == 'model1' | model == 'model2' | model == 'model4') {

      proposedRho <- rhoSampler(model = model,
                                data = data,
                                currentValues = currentValues,
                                proposalSD = proposalSD[[startRho]],
                                currentRho = startRho,
                                priorParameters = priorParameters[[startRho]],
                                nAccepted = nAccepted[startRho])

      # Update the currentValues vector and the nAccepted vector with the new quantities
      currentValues <- proposedRho$currentValues
      nAccepted[startRho] <- proposedRho$nAccepted

    }


    # Fill in the t column of the MCMatrix with the accepted values from the beta, sigma, and rho samplers
    MCMatrix[, t] <- currentValues

  }

  # reduce the size of the MCMatrix according to the burnIn and thinned arguments.
  # In the from argument of the seq function the + 1 is to keep the number of columns in  the MCMatrix the same as the number entered into the thinned argument if a zero is entered into the burnIn argument.
  reducedMCMatrix <- MCMatrix[, seq(from = (nIterations * burnIn + 1), to = nIterations, length.out = thinned)]
  
  
  return (list(MCMatrix = reducedMCMatrix,
               nAccepted = nAccepted))

}
