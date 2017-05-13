# BayesFactorSC is a function that approximates Bayes Factor using the Schwarz Criterion.
# It returns the Bayes Factor approximation in either the log or normal scale.
# The arguments for the function are:
# firstModel - A string indicating which model will be in the numerator of Bayes Factor
# secondModel - A string indicating which model will be in the denominator of Bayes Factor
# firstMHOutput - The first element in the list of the output from the MHSampler function. This is the matrix containing the Markov chain for lambda.
# secondMHOutput - The first element in the list of the output from the MHSampler function. This is the matrix containing the Markov chain for lambda.
# data - A matrix with three columns. The first column contains the gene expression data for T1, the second contains the gene expression data for T2, and the third contains the genetic variant values, V.
# log - A logical argument indicating wether the results are reported on the log scale or the normal scale.
BayesFactorSC <- function (firstModel = 'model0',
                           secondModel = 'model1',
                           firstMHOutput,
                           secondMHOutput,
                           data,
                           log = FALSE) {
  
  # This switch function is used to determine the number of parameters in the model specified by the argument firstModel.
  switch (firstModel,
          
          model0 = {firstD <- 5},
          
          model1 = {firstD <- 6},
          
          model2 = {firstD <- 6},
          
          model3 = {firstD <- 6},
          
          model4 = {firstD <- 7}
          
  )
  
  # This switch function is used to determine the number of parameters in the model specified by the argument secondModel.
  switch (secondModel,
          
          model0 = {secondD <- 5},
          
          model1 = {secondD <- 6},
          
          model2 = {secondD <- 6},
          
          model3 = {secondD <- 6},
          
          model4 = {secondD <- 7}
          
  )
  
  # Calculate the log likelihood for the first model.
  # the first model will be in the numerator of Bayes factor.
  # the arguments that correspond with firstLikelihood are firstModel and firstMHOutput
  firstLikelihood <- calcLogLikelihood(model = firstModel,
                                       data = data,
                                       currentValues = apply(firstMHOutput, 1, median))
  
  # Calculate the log likelihood for the second model.
  # The second model will be in the denominator of Bayes factor.
  # The arguments that correspond with secondLikelihood are secondModel and secondMHOutput
  secondLikelihood <- calcLogLikelihood(model = secondModel,
                                        data = data,
                                        currentValues = apply(secondMHOutput, 1, median))
  
  # Calculate BIC for both the first and second model. length(data[, 1]) gives the sample size of the data.
  firstBIC <- -2 * firstLikelihood + firstD * log(length(data[, 1]))
  secondBIC <- -2 * secondLikelihood + secondD * log(length(data[, 1]))
  
  # Calculate Scharz Criterion using the BIC from the two models
  schwarzCriterion <- -(1 / 2) * (firstBIC - secondBIC)
  
  # Change the scale depending on the input of the argument log
  if (log == TRUE){
    
    return (schwarzCriterion)
    
  } else {
    
    return (exp(schwarzCriterion))
  }
  
}