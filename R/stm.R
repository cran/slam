
## CB 2009/5,6

## simple triple matrix C interfaces
.sums_simple_triplet_matrix <- 
function(x, DIM, na.rm) {
    x <- structure(x$v, i = x$i, j = x$j, Dim = c(x$nrow, x$ncol),
	Dimnames = x$dimnames, class = "stm")
    .Call("_sums_stm", x, DIM, na.rm)

}   

.means_simple_triplet_matrix <-
function(x, DIM, na.rm) {
    x <- structure(x$v, i = x$i, j = x$j, Dim = c(x$nrow, x$ncol),
	Dimnames = x$dimnames, class = "stm")
    s <- .Call("_sums_stm", x, DIM, na.rm)
    if (na.rm) {
	## FIXME inefficient
	x[] <- is.na(x)
	s /(attr(x, "Dim")[-DIM] - .Call("_sums_stm", x, DIM, na.rm))
    } else
	s / attr(x, "Dim")[-DIM]
}


## R interfaces

rowSums <-
function(x, na.rm = FALSE, dims = 1, ...)
   UseMethod("rowSums") 

rowSums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::rowSums(x, na.rm, dims, ...)

rowSums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .sums_simple_triplet_matrix(x, DIM = 1L, na.rm)

colSums <-
function(x, na.rm = FALSE, dims = 1, ...)
   UseMethod("colSums") 

colSums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::colSums(x, na.rm, dims, ...)

colSums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .sums_simple_triplet_matrix(x, DIM = 2L, na.rm)

rowMeans <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("rowMeans")

rowMeans.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::rowMeans(x, na.rm, dims, ...)

rowMeans.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .means_simple_triplet_matrix(x, DIM = 1L, na.rm)

colMeans <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("colMeans")

colMeans.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::colMeans(x, na.rm, dims, ...)

colMeans.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .means_simple_triplet_matrix(x, DIM = 2L, na.rm)


## NOTE the C code must always check for special values and
##      therefore has control over how to proceed. For now
##      it calls the bailout function below.
tcrossprod.simple_triplet_matrix <-
function(x, y = NULL) {
    if (!is(x, "simple_triplet_matrix"))
	stop("'x' not of class simple_triplet_matrix")
    x <- structure(x$v, i = x$i, j = x$j, Dim = c(x$nrow, x$ncol),
	Dimnames = x$dimnames, class = "stm")
    .Call("tcrossprod_stm_matrix", x, y,
	 environment(tcrossprod.simple_triplet_matrix))
}

.tcrossprod.bailout <-
function(x, y) {
    t <- array(0, dim = attr(x, "Dim"), dimnames = attr(x, "Dimnames"))
    t[cbind(attr(x, "i"), attr(x, "j"))] <- x
    tcrossprod(t, y)
}

###
