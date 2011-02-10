
## NOTE this looks all very nice but actually is very
##      inefficient, especially when changing the attributes 
##      of an object the object itself must be copied.
##
## FIXME drop
##

rollup <- 
function(x, MARGIN, INDEX, FUN, ...)
    UseMethod("rollup")

rollup.array <-
function(x, MARGIN, INDEX, FUN = sum, ...) {
    if (is.character(MARGIN))
        MARGIN <- match(MARGIN, names(dimnames(x)))
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
        stop("'MARGIN' invalid")
    if (is.atomic(INDEX))
        INDEX <- list(INDEX)
    if (length(INDEX) != length(MARGIN))
        stop("'INDEX' invalid length")
    names(INDEX) <- MARGIN
    i <- arrayInd(seq_along(x), .dim = dim(x))
    i <- apply(i, 2L, list)
    i <- unlist(i, recursive = FALSE)
    names(i) <- names(dimnames(x))
    for (k in MARGIN) {
	f <- factor(INDEX[[as.character(k)]])
	if (length(f) != dim(x)[k])
	    stop(gettextf("INDEX [%s] invalid length", k))
	i[[k]] <- f[i[[k]]]
    }
    i <- lapply(i, factor)
    x <- lapply(split(x, i), FUN, ...)
    if (all(unlist(lapply(x, length)) == 1L))
	x <- unlist(x, recursive = FALSE)
    i <- lapply(i, levels)
    f <- unlist(lapply(i, length))
    i[-MARGIN] <- list(NULL)
    array(
	data     = x, 
	dim      = f,
	dimnames = i
    )
}

rollup.matrix <- rollup.array

rollup.simple_sparse_array <-
function(x, MARGIN, INDEX, FUN = sum, ...) {
    if (is.character(MARGIN)) 
	MARGIN <- match(MARGIN, names(dimnames(x)))
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
	stop("'MARGIN' invalid")
    if (is.atomic(INDEX))
	INDEX <- list(INDEX)
    if (length(INDEX) != length(MARGIN))
	stop("'INDEX' invalid length")
    names(INDEX) <- MARGIN
    FUN <- match.fun(FUN)
    if (identical(FUN, sum)) {
	for (k in MARGIN) {
	    i <- factor(INDEX[[as.character(k)]])
	    if (length(i) != dim(x)[k])
		stop(gettextf("INDEX [%s] invalid length", k))
	    x$dim[k]       <- length(levels(i))
	    dimnames(x)[k] <- list(levels(i))
	    i <- c(i)[x$i[, k, drop = TRUE]]
	    x$i[, k] <- i 
	    i <- is.na(i)
	    if (any(i)) {
		i <- !i
		x$i <- x$i[i, ,drop = FALSE]
		x$v <- x$v[i]
	    }
	    i <- apply(x$i, 1L, paste, collapse = "\r")
	    x$v <- lapply(split(x$v, i), FUN, ...)
	    x$i <- x$i[match(names(x$v), i),, drop = FALSE]
	    x$v <- unlist(x$v, use.names = FALSE)
	}
	x
    } else 
	stop("'FUN' not implemented")
}

rollup.simple_triplet_matrix <- 
function(x, MARGIN, INDEX, FUN = sum, ...) {
    FUN <- match.fun(FUN)
    if (identical(FUN, sum)) {
	if (is.character(MARGIN)) 
	    MARGIN <- match(MARGIN, names(dimnames(x)))
	if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
	    stop("'MARGIN' invalid")
	if (is.atomic(INDEX))
	    INDEX <- list(INDEX)
	if (length(INDEX) != length(MARGIN))
	    stop("'INDEX' invalid length")
	names(INDEX) <- MARGIN
	for (k in MARGIN)
	    x <- switch(k,
		{
		    x <- t(x)
		    x <- rollup(x, 2L, INDEX[as.character(k)], FUN, ...)
		    t(x)
		},
		.Call("_row_tsums", 
		      x, factor(INDEX[[as.character(k)]]), 
		      if (is.null(list(...)$na.rm))
			  FALSE
		      else
			  as.logical(list(...)$na.rm), 
		      FALSE),
		stop("'MARGIN' invalid")
	    )
	x
    } else
	stop("'FUN' not implemented")
}

##
rollup.default <-
function(x, MARGIN, INDEX, FUN = sum, ...) {
    if (length(dim(x)))
	stop("dim(x) must have a positive length")
    x <- as.array(x)
    rollup(x, MARGIN, INDEX, FUN, ...)
}

###
