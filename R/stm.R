
## CB 2009/5,6,10 2010/6


.means_simple_triplet_matrix <-
function(x, DIM, na.rm)
{
    s <- .Call("_sums_stm", x, DIM, na.rm)
    n <- c(x$nrow, x$ncol)[-DIM]
    if (na.rm) {
	x$v <- is.na(x$v)
	nna <- .Call("_sums_stm", x, DIM, FALSE)
	s / (n - nna)
    }
    else
	s /  n
}


## R interfaces

row_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("row_sums")

row_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::rowSums(x, na.rm, dims, ...)

row_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .Call("_sums_stm", x, 1L, na.rm)

col_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("col_sums")

col_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::colSums(x, na.rm, dims, ...)

col_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .Call("_sums_stm", x, 2L, na.rm)

row_means <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("row_means")

row_means.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::rowMeans(x, na.rm, dims, ...)

row_means.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .means_simple_triplet_matrix(x, DIM = 1L, na.rm)

col_means <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("col_means")

col_means.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base:::colMeans(x, na.rm, dims, ...)

col_means.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .means_simple_triplet_matrix(x, DIM = 2L, na.rm)


## NOTE the C code must always check for special values and
##      therefore has control over how to proceed. For now
##      it calls the bailout function below. For verbose
##      information set the last argument to TRUE.
##
##      The symmetric case is now also handled in C. Runtime
##      could be further improved if data need not to be 
##      ordered (see the C code).
tcrossprod_simple_triplet_matrix <-
function(x, y = NULL) {
    if (is.null(y))
	.Call("tcrossprod_stm_stm", x, y, 
	      environment(tcrossprod_simple_triplet_matrix), FALSE)
    else
	.Call("tcrossprod_stm_matrix", x, y,
	      environment(tcrossprod_simple_triplet_matrix), FALSE, FALSE)
}

## For now internal.
.ttcrossprod_simple_triplet_matrix <-
function(x, y = NULL) {
    if (is.null(y))
	tcrossprod_simple_triplet_matrix(x)
    else
	.Call("tcrossprod_stm_matrix", x, y,
	      environment(tcrossprod_simple_triplet_matrix), FALSE, TRUE)
}

## FIXME warning?
.tcrossprod.bailout <-
function(x, y, transpose) {
    if (transpose) {
	if (is.null(y))
	    tcrossprod(as.matrix(x))
	else
	    tcrossprod(y, as.matrix(x))
    }
    else
	tcrossprod(as.matrix(x), y)
}

##
.nnzero <- 
function(x, scale = FALSE) {
    v <- c("simple_triplet_matrix", "simple_sparse_array")
    if (inherits(x, v))
	v <- x$v
    else {
	x <- as.array(x)
	v <- x
    }
    v <- v == vector(typeof(v), 1L)
    v <- v + 1L
    n <- length(v)
    v <- tabulate(v, 2L)
    v <- c(v, n - sum(v))
    names(v) <- c("nnzero", "nzero", NA)
    if (scale)
	v <- v / prod(dim(x))
    v
}

###
