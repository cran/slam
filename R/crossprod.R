
## NOTE the C code must always check for special values and
##      therefore has control over how to proceed. For now
##      it calls the bailout function below. 
##
##	For verbose information set the verbose argument to 
##	TRUE. Transposition of the return value (!) is only
##	implemented for dense.
##
##      The general case is now also handled in C. Runtime
##      could be further improved if the data need not be 
##      ordered (see the C code).

.tcrossprod_simple_triplet_matrix <-
function(x, y = NULL, transpose = FALSE, bailout = TRUE, verbose = FALSE) {
    if (!is.simple_triplet_matrix(x))
	stop("'x' not of class simple_triplet_matrix")
    if (is.null(y) ||
	is.simple_triplet_matrix(y)) {
	if (transpose)
	    stop("'transpose' not implemented")
	.Call(R_tcrossprod_stm_stm, x, y,
	    if (bailout)
		environment(.tcrossprod_simple_triplet_matrix),
	    verbose
	)
    }
    else
	.Call(R_tcrossprod_stm_matrix, x,
	    as.matrix(y),
	    if (bailout)
		environment(.tcrossprod_simple_triplet_matrix),
	    verbose,
	    transpose
	)

}

.tcrossprod_bailout <-
function(x, y, transpose) {
    if (transpose)
	## see above
	base::tcrossprod(y, as.matrix(x))
    else
	base::tcrossprod(as.matrix(x),
	    if (is.null(y))
		y
	    else
		as.matrix(y)
	)
}

## Used by package skmeans.
.ttcrossprod_simple_triplet_matrix <-
function(x, y = NULL)
    .tcrossprod_simple_triplet_matrix(x, y, TRUE)

##
tcrossprod_simple_triplet_matrix <-
function(x, y = NULL) {
    if(is.simple_triplet_matrix(x)) {
        if(!is.simple_triplet_matrix(y) &&
           .is_sparse_mat_coercible_to_stm(y))
            y <- as.simple_triplet_matrix(y)
        .tcrossprod_simple_triplet_matrix(x, y)
    }
    else if(is.simple_triplet_matrix(y)) {
        x <- if(.is_sparse_mat_coercible_to_stm(x))
                 as.simple_triplet_matrix(x)
             else
                 as.matrix(x)
        .tcrossprod_simple_triplet_matrix(y, x, TRUE)
    }
    else
        stop("neither 'x' nor 'y' of class 'simple_triplet_matrix'")
}

crossprod_simple_triplet_matrix <-
function(x, y = NULL) {
    if(is.simple_triplet_matrix(x)) {
        y <- if(is.null(y))
                 y
             else if(is.simple_triplet_matrix(y))
                 t(y)
             else if(.is_sparse_mat_coercible_to_stm(y))
                 t(as.simple_triplet_matrix(y))
             else
                 t(as.matrix(y))
        .tcrossprod_simple_triplet_matrix(t(x), y)
    }
    else if(is.simple_triplet_matrix(y)) {
        x <- if(.is_sparse_mat_coercible_to_stm(x))
                 as.simple_triplet_matrix(x)
             else
                 as.matrix(x)
        .tcrossprod_simple_triplet_matrix(t(y), t(x), TRUE)
    }
    else
        stop("neither 'x' nor 'y' of class 'simple_triplet_matrix'")
}

matprod_simple_triplet_matrix <-
function(x, y) {
    if(is.simple_triplet_matrix(x)) {
        y <- if(is.simple_triplet_matrix(y))
                 y
             else if(.is_sparse_mat_coercible_to_stm(y))
                 as.simple_triplet_matrix(y)
             else
                 as.matrix(y)
	.tcrossprod_simple_triplet_matrix(x, t(y))
    }
    else if(is.simple_triplet_matrix(y)) {
        x <- if(.is_sparse_mat_coercible_to_stm(x))
                 as.simple_triplet_matrix(x)
             else
                 as.matrix(x)
        .tcrossprod_simple_triplet_matrix(t(y), x, TRUE)
    }
    else
        stop("neither 'x' nor 'y' of class 'simple_triplet_matrix'")
}

##

matrixOps.simple_triplet_matrix <-
function(x, y)
{
    switch(.Generic,
           "%*%" = matprod_simple_triplet_matrix(x, y),
           "crossprod" = if(missing(y))
                             crossprod_simple_triplet_matrix(x)
                         else
                             crossprod_simple_triplet_matrix(x, y),
           "tcrossprod" = if(missing(y))
                              tcrossprod_simple_triplet_matrix(x)
                          else
                              tcrossprod_simple_triplet_matrix(x, y))
}

chooseOpsMethod.simple_triplet_matrix <- 
function(x, y, mx, my, cl, reverse)
    TRUE
