\name{MHSampler}
\alias{MHSampler}

\title{
Metropolis Hastings Algorithm
}

\description{
Runs a Metropolis Hastings algorithm on a given set of data
}

\usage{
MHSampler(nIterations = 100, model = 'model0',
          burnIn = 0.2,
          thinned = 300,
          data,
          parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72),
          priorParameters = list(b0.1 = c(0, 1), b1.1 = c(0, 1), b0.2 = c(0, 1), b1.2 = c(0, 1), sd.1 = c(1, 1), sd.2 = c(1, 1), rho = c(0, 0.5)),
          proposalSD = list(b0.1 = 1, b1.1 = 1, b0.2 = 1, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.5))
}

\arguments{
  \item{nIterations}{
  Positive integer. The number of times to run the sampler
  }

  \item{burnIn}{
  Real number between 0 and 1. The percentage of the Markov chain that will be discarded.
  }

  \item{thinned}{
  Positive integer. The number that the Markov chain will be thinned to.
  }

  \item{model}{
  There are five different models. To select a certain model use the form 'modelX' where X is an integer between 0 and 4. The default is 'model0'.
  }

  \item{data}{
  Data frame. First column must be the gene expression data for gene one, the second column must be the gene expression data for gene two, and the third column must be the genotype data.
  }

  \item{parameters}{
  List of seven. This list is used as the starting values for the sampler. The elements of the list are:

  b0.1 -
    Intercept in the linear model b0.1 + b1.1 * V. This linear model is the mean of the expression of gene 1.

    b1.1 -
    Slope in the linear model b0.1 + b1.1 * V. This linear model is the mean of the expression of gene 1.

    b0.2 -
    Intercept in the linear model b0.2 + b1.2 * V. This linear model is the mean of the expression of gene 2.

    b1.2 -
    Slope in the linear model b0.2 + b1.2 * V. This linear model is the mean of the expression of gene 2.

    sd.1 -
    Standard deviation for the gene T1.

    sd.2 -
    Standard deviation for the gene T2.

    rho -
    Correlation between the two genes.
  }

  \item{priorParameters}{
  List of seven. This is a list of means and standard deviations for the prior distribution of each parameter.
  }

  \item{proposalSD}{
  List of seven. This list is the standard deviation for the proposal distribution. All proposal distributions are normal. In the case of the standard deviation and correlation parameters a truncated normal distribution is used.
  }
}

\value{
Returns a list of two. The first element is the Markov Chain for each parameter. The second element is a vector that shows how many times the proposed value for each parameter was accepted.
}

\examples{
# Values used to simulate the data for model 3.
actualValues <- list(b0.1 = 3.2, b1.1 = 1.7, b0.2 = 4.6, b1.2 = 3.1, sd.1 = 2.4, sd.2 = 0.9, rho = 0.72)

# Number of iterations to run the MHSampler function.
m <- 5000

Model3 <- simulateData(N = 100,
                       p = 0.43,
                       seed = 338,
                       model = 'model3',
                       parameters = actualValues)

Output3 <- MHSampler(nIterations = m,
                     burnIn = 0.2,
                     thinned = 200,
                     model = 'model0',
                     data = Model3,
                     parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72),
                     priorParameters = list(b0.1 = c(0, 1), b1.1 = c(0, 1), b0.2 = c(0, 1), b1.2 = c(0, 1), sd.1 = c(1, 1), sd.2 = c(1, 1), rho = c(0, 0.5)),
                     proposalSD = list(b0.1 = 1, b1.1 = 1, b0.2 = 1, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.5))
}
