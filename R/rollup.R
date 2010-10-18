
## NOTE this looks all very nice but actually is very
##      inefficient especially as changing the attributes 
##      of an object the object itself must be copied.
##
## TBDS 1) how can we reliably test if FUN is scalar, e.g.
##	   length(FUN(1:2)) == 1L?
##      2) should we map character, e.g.
##         if (is.character(MARGIN))
##	       MARGIN <- match(MARGIN, names(dimnames(x)))

rollup <- 
function(x, MARGIN, INDEX, FUN, ...)
    UseMethod("rollup")

rollup.array <-
function(x, MARGIN, INDEX, FUN, ...) {
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
	stop("'MARGIN' invalid")
    if (is.atomic(INDEX))
	INDEX <- list(INDEX)
    if (length(INDEX) != length(MARGIN))
	stop("'INDEX' invalid length")
    names(INDEX) <- MARGIN
    for (k in MARGIN) {
	m <- seq_along(dim(x))[-k]
	f <- factor(INDEX[[as.character(k)]])
	d <- dimnames(x)
	d[k] <- list(levels(f))
	x <- array(
	    apply(x, m, tapply, f, FUN, ...),
	    dim      = c(length(levels(f)), dim(x)[m]), 
	    dimnames = d[c(k, m)]
	)
	x <- aperm(x, perm = order(c(k, m)))
    }
    x
}

rollup.matrix <- rollup.array

rollup.simple_sparse_array <-
function(x, MARGIN, INDEX, FUN, ...) {
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
	stop("'MARGIN' invalid")
    if (is.atomic(INDEX))
	INDEX <- list(INDEX)
    if (length(INDEX) != length(MARGIN))
	stop("'INDEX' invalid length")
    names(INDEX) <- MARGIN
    for (k in MARGIN) {
	i <- factor(INDEX[[as.character(k)]])
	if (length(i) != dim(x)[k])
	    stop(gettextf("INDEX [%s] invalid length", k))
	x$i[, k] <- l <- c(i)[x$i[, k, drop = TRUE]]
	l <- !is.na(l)
	if (!all(l)) {
	    x$i <- x$i[l, ,drop = FALSE]
	    x$v <- x$v[l]
	}
	ind <- apply(x$i, 1L, paste, collapse = ".")
	x$v <- c(tapply(x$v, ind, FUN, ...))
	x$i <- x$i[match(names(x$v), ind),, drop = FALSE]
	x$dim[k]         <- length(levels(i))
	dimnames(x)[k]   <- list(levels(i))
    }
    names(x$v) <- NULL
    x
}

rollup.simple_triplet_matrix <- 
function(x, MARGIN, INDEX, FUN, ...) {
    FUN <- match.fun(FUN)
    if (identical(FUN, sum)) {
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
    } 
    else {
	x <- as.simple_sparse_array(x)
	x <- rollup(x, MARGIN, INDEX, FUN, ...)
	as.simple_triplet_matrix(x)
    }
}

##
rollup.default <-
function(x, MARGIN, INDEX, FUN, ...) {
    if (length(dim(x)))
	stop("dim(x) must have a positive length")
    x <- as.array(x)
    rollup(x, MARGIN, INDEX, FUN, ...)
}

###
