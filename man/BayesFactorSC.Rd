\name{BayesFactorSC}
\alias{BayesFactorSC}

\title{
Bayes Factor Approximation using Schwarz Criterion
}

\description{
Approximates Bayes Factor using Schwarz Criterion. The function uses BIC to calculate Schwarz Criterion.
}

\usage{
BayesFactorSC(firstModel = 'model0',
              secondModel = 'model1',
              firstMHOutput,
              secondMHOutput,
              data,
              log = FALSE)
}

\arguments{
  \item{firstModel}{
  String indicating the model that will be in the numerator of Bayes Factor.
  }

  \item{secondModel}{
  String indicating the model that will be in the denominator of Bayes Factor.
  }

  \item{firstMHOutput}{
  Markov chain matrix of the parameters from the model specified in the firstModel argument. This is the first element in the list of the output from the MHSampler function.
  }

  \item{secondModel}{
  Markov chain matrix of the parameters from the model specified in the secondModel argument. This is the first element in the list of the output from the MHSampler function.
  }

  \item{data}{
  Data frame. First column must be the gene expression data for gene one, the second column must be the gene expression data for gene two, and the third column must be the genotype data.
  }

  \item{log}{
  Logical argument. If TRUE the Bayes Factor approximation will be returned on the log scale.
  }
}

\value{
A real number. The scale depends on the input of the log argument.
}

\references{
Robert E. Kass and Adrian E. Raftery.
\emph{Bayes Factors}.
Journal of the American Statistical Association, Vol.90, 773-795, 1995.
}

\examples{
# Values used to simulate the data for model 3.
actualValues <- list(b0.1 = 3.2, b1.1 = 1.7, b0.2 = 4.6, b1.2 = 3.1, sd.1 = 2.4, sd.2 = 0.9, rho = 0.72)

# Number of iterations to run the MHSampler function.
m <- 5000

# Simulate data under model 3
Model3 <- simulateData(N = 100,
                       p = 0.43,
                       seed = 338,
                       model = 'model3',
                       parameters = actualValues)

# Run the MHSampler function for models 3 and 4 and calculate Bayes Factor for the two models.
Output3 <- MHSampler(nIterations = m,
                     burnIn = 0.2,
                     thinned = 200,
                     model = 'model3',
                     data = Model3,
                     parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72),
                     priorParameters = list(b0.1 = c(0, 1), b1.1 = c(0, 1), b0.2 = c(0, 1), b1.2 = c(0, 1), sd.1 = c(1, 1), sd.2 = c(1, 1), rho = c(0, 0.5)),
                     proposalSD = list(b0.1 = 1, b1.1 = 1, b0.2 = 1, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.5))

Output4 <- MHSampler(nIterations = m,
                     burnIn = 0.2,
                     thinned = 200,
                     model = 'model4',
                     data = Model3,
                     parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72),
                     priorParameters = list(b0.1 = c(0, 1), b1.1 = c(0, 1), b0.2 = c(0, 1), b1.2 = c(0, 1), sd.1 = c(1, 1), sd.2 = c(1, 1), rho = c(0, 0.5)),
                     proposalSD = list(b0.1 = 1, b1.1 = 1, b0.2 = 1, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.5))

BFA <- BayesFactorSC(firstModel = 'model3',
                     secondModel = 'model4',
                     firstMHOutput = Output3[[1]],
                     secondMHOutput = Output4[[1]],
                     data = Model3,
                     log = FALSE)
}
