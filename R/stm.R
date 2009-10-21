
## CB 2009/5,6,10


.as_stm.simple_triplet_matrix <- 
function(x, compat = FALSE) {
    if (!inherits(x, "simple_triplet_matrix"))
	stop("'x' not of class simple_triplet_matrix")
    ## old-style
    if(compat)
	.structure(x$v, i = x$i, j = x$j, Dim = c(x$nrow, x$ncol),
                   Dimnames = x$dimnames, class = "stm")
    else
	.Call("_as_stm", x$v, x$i, x$j, c(x$nrow, x$ncol), x$dimnames)
}

## simple triple matrix C interfaces
.sums_simple_triplet_matrix <- 
function(x, DIM, na.rm)
{
    x <- .as_stm.simple_triplet_matrix(x)
    .Call("_sums_stm", x, DIM, na.rm)
}   

.means_simple_triplet_matrix <-
function(x, DIM, na.rm)
{
    x <- .as_stm.simple_triplet_matrix(x)
    s <- .Call("_sums_stm", x, DIM, na.rm)
    if (na.rm) {
	## FIXME inefficient
	if (is.list(x))
	    attr(x, "v")[] <- is.na(attr(x, "v"))
	else
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
function(x, y = NULL)
{
    x <- .as_stm.simple_triplet_matrix(x)
    .Call("tcrossprod_stm_matrix", x, y,
          environment(tcrossprod.simple_triplet_matrix))
}

.tcrossprod.bailout <-
function(x, y)
{
    t <- array(0, dim = attr(x, "Dim"), dimnames = attr(x, "Dimnames"))
    t[cbind(attr(x, "i"), attr(x, "j"))] <- 
	if(is.list(x)) attr(x, "v") else x
    tcrossprod(t, y)
}

## KH 2009/10

.structure <-
function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))

###
