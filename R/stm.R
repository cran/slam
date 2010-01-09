
## CB 2009/5,6,10


.means_simple_triplet_matrix <-
function(x, DIM, na.rm)
{
    s <- .Call("_sums_stm", x, DIM, na.rm)
    if (na.rm) {
	x$v <- is.na(x$v)
	s /(c(x$nrow, x$ncol)[-DIM] - .Call("_sums_stm", x, DIM, na.rm))
    } else
	s / c(x$nrow, x$ncol)[-DIM]
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
tcrossprod_simple_triplet_matrix <-
function(x, y = NULL)
    .Call("tcrossprod_stm_matrix", x, y,
          environment(tcrossprod_simple_triplet_matrix), FALSE)

## FIXME warning?
.tcrossprod.bailout <-
function(x, y)
    tcrossprod(as.matrix(x), y)

###
