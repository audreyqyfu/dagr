# correlation calculates the correlation between two vectors of gene expression levels.
# It returns a dataframe with four columns. The first column has the correlation between the two gene expression level vectors. The following three columns contain the correlation between the two gene expression level vectors when the value of the genotype vector is 0, 1, and 2.
# The arguments for the function are:
# t1 - A vector of gene expression levels.
# t2 - A vector of gene expression levels.
# V - A vector of integers representing the genotype of the individual.
correlation <- function (t1, t2, V) {

  cAll <- cor(t1, t2)
  cors <- list()

  for (i in 1:3) {

    cors[[i]] <- cor(t1[which(V == (i - 1))],
                     t2[which(V == (i - 1))])

  }

  return (data.frame(CorAll = cAll,
                     Cor0 = cors[[1]],
                     Cor1 = cors[[2]],
                     Cor2 = cors[[3]]))
}
