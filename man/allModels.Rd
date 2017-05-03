\name{allModels}
\alias{allModels}

\title{
All Models
}

\description{
Simulates expression levels for all five models at once.
}

\usage{
allModels(N, p, seed, parameters)
}

\arguments{
  \item{N}{
  Positive integer. The number of observations to simulate.
  }

  \item{p}{
  Real number between 0 and 1. Proportion of the population that has the reference allele.
  }

  \item{seed}{
  A positive number used for random number generation.
  }

  \item{parameters}{
  A list containing seven elements:

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
}

\value{
This function returns a dataframe with eleven columns. The first ten columns contain the data for the expression for genes 1 and 2, for the five different models. The 11th column contains the genotype of the eQTL. The data frame that is returned will have some columns with identical output. This is because the seed is set at the beginning of the function and some topologies share the same interaction between the eQTL and the two genes.
}

\examples{
data <- allModels(N = 100,
                  p = 0.43,
                  seed = 338,
                  parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72))
}
