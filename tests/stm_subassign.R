
##
library("slam")

s <- as.simple_triplet_matrix(diag(4))
s[1:8] <- 1:8
as.matrix(s)

s[2:3,] <- 1:8
as.matrix(s)

s[,2:3] <- 1:8
as.matrix(s)

s[] <- 1:8
as.matrix(s)

###
