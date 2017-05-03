\name{Output4}
\alias{Output4}
\docType{data}

\title{
Output from the Metropolis Hastings algorithm
}

\description{
This data set is a list of two. The first element in the list is the Markov Chain for each parameter in model 4. The second element in the list is a vector of lenght seven which contains the number of times the proposed value was accepted for each parameter in model 4.
}

\usage{data("Output4")}

\format{
  A list of two. The first element in the list is a matrix of seven rows and 5000 columns. The rows are the parameters in model four. The columns contain the Markov Chain from the Metropolis Hastings algorithm.
  \describe{
    \item{\code{beta01}}{a numeric vector}
    \item{\code{beta11}}{a numeric vector}
    \item{\code{beta02}}{a numeric vector}
  }
}

\examples{
data(Output4)
## maybe str(Output4) ; plot(Output4) ...
}

\keyword{datasets}
