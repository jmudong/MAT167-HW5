library(MASS)

# Creating hyperlink matrix
G = t(fractions(matrix(c(0, 0, 0, 1/2, 0, 0, 0, 1/2, 0, 0, 0, 0,
                       1/3, 0, 0, 1/3, 1/3, 0, 0, 0, 0, 0, 0, 0,
                       0, 1/2, 0, 0, 1/2, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 1/2, 0, 1/2, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 1/4, 1/4, 1/4, 0, 1/4, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 1/2, 0, 1/2, 0,
                       0, 0, 0, 0, 0, 0, 1/2, 0, 0, 0, 1/2, 0,
                       0, 0, 0, 0, 0, 1/2, 1/2, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 1/3, 0, 0, 0, 1/3, 0, 1/3,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                       nrow = 12, ncol = 12)))

# Eigenvalue Decomposition
ev1 = eigen(t(G))
evalue1 = ev1$values[1] # Choosing largest eigenvalue
evector1 = ev1$vectors[,1] # Eigenvector corresponding to largest eigenvalue

# Adding artificial links to dangling node
Gnew = G
Gnew[5,] = rep(1/12, 12)
ev2 = eigen(t(Gnew)) # Eigenvalue decomposition
evalue2 = ev2$values[1]
evector2 = ev2$vectors[,1]


