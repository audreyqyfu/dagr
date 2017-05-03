\name{correlation}
\alias{correlation}

\title{
Correlation
}

\description{
Calculates the correlation between two gene expression vectors at all levels of the genotype. In addition it calculates the correlation between the two genes at each level of the genotype.
}

\usage{
correlation(t1, t2, V)
}

\arguments{
  \item{t1}{
  Gene expression vector.
  }

  \item{t2}{
  Gene expression vector.
  }

  \item{V}{
  Genotype vector.
  }
}

\value{
This function returns a dataframe where the first column is the correlation of T1 and T2 for all values of V, the second column is the correlation of T1 and T2 when V = 0, the third column is the correlation of T1 and T2 when V = 1, and the fourth column is the correlation of T1 and T2 when V = 2.
}

\examples{
Model1 <- simulateData(N = 100,
                       p = 0.43,
                       seed = 338,
                       model = 'model1',
                       parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72))

correlation(Model1[, 1], Model1[, 2], Model1[, 3])
}
