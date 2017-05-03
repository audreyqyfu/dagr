\name{scatterPlot}
\alias{scatterPlot}

\title{
Scatter Plot
}

\description{
Creates a scatter plot of two gene expression vectors. It also displays the mean of the two gene expression vectors at each level of the genotype.
}

\usage{
scatterPlot(t1, t2, V)
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
Creates a scatter plot of the data.
}

\examples{
Model1 <- simulateData(N = 100,
                       p = 0.43,
                       seed = 338,
                       model = 'model1',
                       parameters = list(b0.1 = 0, b1.1 = 1, b0.2 = 0, b1.2 = 1, sd.1 = 1, sd.2 = 1, rho = 0.72))

scatterPlot(Model1[, 1], Model1[, 2], Model1[, 3])
}
