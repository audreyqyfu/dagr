\name{MHPlots}
\alias{MHPlots}

\title{
Create MCMC Plots
}

\description{
Creates a line plot of the Markov Chain for each parameter. It also creates a histogram of the Markov Chain with a vertical line of the actual value used to generate the data for each parameter.
}

\usage{
MHPlots(iterations, model, MHOutput, actualValues, parameter, first, last)
}

\arguments{

  \item{iterations}{
  Positive interger. This is the number of times the MHSampler function was run.
  }

  \item{model}{
  There are five different models. To select a certain model use the form 'modelX' where X is an integer between 0 and 4. The default is 'model0'.
  }

  \item{MHOutput}{
  This is the output from the MHSampler function. It is a list containing the matrix of Markov Chains for each parameter and the number of times the proposed value for each parameter was selected.
  }

  \item{actualValues}{
  A list containing the values for each of the parameters that was used to simulate the data. The elements of the list are:

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

  \item{parameter}{
  An interger representing the location in the list of parameters for which the plot will be created.
  }

  \item{first}{
  Positive integer. Location in the Markov Chain where the values will start to be used for the line plot and the histogram. Usually this value will be the number of times the Metropolis Hastings sampler was run multiplied by 0.2.
  }

  \item{last}{
  Positive integer. Location in the Markov Chain to stop taking values to include in the line plot and the histogram. This value should be the number of times the Metropolis Hastings sampler was run.
  }

}

\value{
This function returns a line plot and a histogram for each parameter.
}

\examples{
data(Output4)

# Values used to simulate the data for model 4.
actualValues <- list(b0.1 = 3.2, b1.1 = 1.7, b0.2 = 4.6, b1.2 = 3.1, sd.1 = 2.4, sd.2 = 0.9, rho = 0.72)

# Number of times to run the MHSampler function.
m <- 5000

# Create the plots for the parameters in model 4
for (i in 1:length(actualValues)) {

  MHPlots(iterations = m,
          model = 'model4',
          MHOutput = Output4,
          actualValues = actualValues,
          parameter = i,
          first = m * 0.2,
          last = m)

}
}
