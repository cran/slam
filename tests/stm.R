
library("slam")
set.seed(20090626)

###

x <- sample(0:5, 100, T, prob=c(.8,rep(.04,5)))
x <- matrix(as.logical(x), nrow = 20,
     dimnames = list(rows = 1:20, cols = LETTERS[1:5]))
x

xst <- as.simple_triplet_matrix(x)
xst

identical(rowSums(x), rowSums(xst))
identical(colSums(x), colSums(xst))
identical(rowMeans(x), rowMeans(xst))
identical(colMeans(x), colMeans(xst))

## NAs

xna <- x
n <- prod(dim(x))
is.na(xna) <- sample(seq_len(n), ceiling(n * .1))
xna

xnast <- as.simple_triplet_matrix(xna)
xnast

identical(rowSums(xna), rowSums(xnast))
identical(colSums(xna), colSums(xnast))
identical(rowMeans(xna), rowMeans(xnast))
identical(colMeans(xna), colMeans(xnast))

identical(rowSums(xna, na.rm = TRUE), rowSums(xnast, na.rm = TRUE))
identical(colSums(xna, na.rm = TRUE), colSums(xnast, na.rm = TRUE))
identical(rowMeans(xna, na.rm = TRUE), rowMeans(xnast, na.rm = TRUE))
identical(colMeans(xna, na.rm = TRUE), colMeans(xnast, na.rm = TRUE))

## cross-product

identical(tcrossprod(x), tcrossprod.simple_triplet_matrix(xst))
identical(tcrossprod(x), tcrossprod.simple_triplet_matrix(xst, x))

x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 3)
x
s <- as.simple_triplet_matrix(x)

identical(tcrossprod(x), tcrossprod.simple_triplet_matrix(s))
identical(tcrossprod(x), tcrossprod.simple_triplet_matrix(s, x))

###


