# scatterPlot creates a scatter plot from two gene expression level vectors.
# The arguments for the function are:
# t1 - A vector of gene expression levels.
# t2 - A vector of gene expression levels.
# V - A vector of integers representing the genotype of the individual.
scatterPlot <- function (t1, t2, V) {

  sPlot <- plot(t1,
                t2,
                col = ifelse(V == 0, 'orange', ifelse(V == 1, 'lightseagreen', 'purple')),
                pch = ifelse(V == 0, 2, ifelse(V == 1, 3, 0)))

  legend(x = 'bottomright',
         legend = c('v = 0', 'v = 1', 'v = 2'),
         pch = c(2, 3, 0),
         col = c('orange', 'lightseagreen','purple'))

  for (e in 1:3) {

    points(mean(t1[which(V == e - 1)]),
           mean(t2[which(V == e - 1)]),
           col = 'black',
           pch = 15)

  }
}
