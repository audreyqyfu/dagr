\name{simulateData}
\alias{simulateData}

\title{
Simulates data for the specified model.
}

\description{
Simulates expression levels for two genes, denoted T1 and T2, and the genotype of the eQTL, denoted V, for a specified model.
}

\usage{
simulateData(N, p, seed = 338, model = 'model0', parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0))
}

\arguments{

  \item{N}{
  A positive scalar for the number of observations to generate.
  }

  \item{p}{
  Population frequency of the reference allele. A real number between 0 and 1.
  }

  \item{seed}{
  A positive number used for random number generation.
  }

  \item{model}{
  There are five different models. To select a certain model use the form 'modelX' where X is an integer between 0 and 4. The default is 'model0'.
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
This function returns a dataframe with three columns. The first two columns contain the expression data for genes 1 and 2. The third column are the data for the genotype.
}

\examples{
Model1 <- simulateData(N = 100,
                       p = 0.43,
                       seed = 338,
                       model = 'model1',
                       parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72))
}
