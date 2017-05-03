# MHPlots creates a line plot and a histogram of the sample created from the MHSampler function
# It returns two plots. The first is the line plot and the second is the histogram
# The arguments of the function are:
# iterations - A scalar indicating the number of times the MHSampler function ran.
# MHOutput - A list of two. The first element of the list is a vector. Each element of the vector represents how many proposed values were accepted for each parameter. The second element of the list is a matrix containing the Markov Chain for each parameter.
# actualValues - A list of length seven. Each element of the list is the parameter value used to generate the data.
# first - An integer indicating at what column to start using for the plots.
# last - An integer indicating at whe column to stop using for the plots.
MHPlots <- function (iterations,
                     model,
                     MHOutput,
                     actualValues,
                     parameter,
                     first,
                     last) {

  plot(MHOutput[[1]][parameter, first:last],
       type = 'l',
       main = c(model, names(actualValues[parameter])),
       xlab = 'Last 80% of Iterations',
       ylab = 'Parameter Value')

  hist(MHOutput[[1]][parameter, first:last],
       main = c(model, names(actualValues[parameter])),
       xlab = 'Parameter Value')
  abline(v = actualValues[[parameter]], lwd = 2, col = 'orange')

}
