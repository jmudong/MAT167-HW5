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
evector2 = abs(ev2$vectors[,1])

# Normalizing eigenvector
evector2n = evector2/sum(evector2)

# Ordering eigenvector
ranks = data.frame(webpage = letters[1:12], pagerank = evector2n)
ranks_sorted = ranks[order(ranks$pagerank, decreasing = TRUE),]

# Power Iteration
powermethod = function(A, k){
  # function that returns largest eigenvalue and eigenvector of matrix A
  # given k iterations
  v = fractions(rep(1/nrow(A), nrow(A))) # initialize vector
  for(i in 1:k){
    w = A %*% v
    v = w/sum(abs(w))
    lambda = t(v) %*% A %*% v
  }
  return(list(eigenvalue = lambda, eigenvector = v))
}

# Iterating through different k's



# Term-Document Matrix
D = t(matrix(c(1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0,
               1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
               1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0,
               1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1,
               0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0,
               0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0,
               1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
               1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
               0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1,
               0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1,
               1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
               1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1,
               0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0,
               0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0,
               0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
               0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1),
             ncol = 16, nrow = 12))
