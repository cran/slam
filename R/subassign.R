## CB 2012/9
##
## FIXME extending might be useful unless implemented
##	 as for dense arrays.
##
`[<-.simple_sparse_array` <- 
function(x, ..., value) {

    if (!length(value))
	stop("replacement has length zero")

    nd <- length(x$dim)
    pd <- prod(x$dim)

    ## Disable features that may exhaust resources.
    .protect <- pd > 16777216L

    na <- nargs()
    if (na == 3L && missing(..1))
	if (.protect)
	    stop("Empty subscripting disabled.")
	else
	    return(
		do.call("[<-.simple_sparse_array",
			list(x = x, seq_len(pd), value = value))
	    )
	    
    ## Single index subscripting.
    if (na == 3L) {
	I <- ..1
	if (!is.numeric(I))
	    stop("Only numeric subscripting is implemented.")
	if (!length(I))
	    return(x)
	## Missing values in subscripts.
	k <- is.na(I)
	if (any(k))
	    if (length(value) == 1L)
		I[k] <- 0L
	    else
		stop("NAs are not allowed in subscripted assignments")
	rm(k)
	## Vector subscripting.
	if (!is.matrix(I)) {
	    ## 52-bit safe
	    if (pd > 4503599627370496)
		stop("Vector subscripting disabled for this object.")
	    ## Map.
	    if (is.double(I))
		I <- trunc(I)
	    if (all(I >= 0L)) {
		## Remove zero subscripts.
		I <- I[I > 0L]
		if (!length(I))
		    return(x)
		if (any(I > pd))
		    stop("Extending is not implemented.")
	    } else {
		if (.protect)
		    stop("Negative subscripting disabled for this object.")
		if (all(I <= 0L)) {
		    ## NOTE this fails if NAs are introduced by 
		    ##	    coercion. 
		    I <- seq_len(pd)[I]
		} else
		    stop("only 0's may be mixed with negative subscripts")
	    }
	    ## Expand.
	    I <- arrayInd(I, .dim = x$dim)
	} else
	    ## NOTE as the other replacement rules are no less 
	    ##	    confusing we allow this, too.
	    if (ncol(I) != nd) {
		dim(I) <- NULL
		return(
		    do.call("[<-.simple_sparse_array", 
			    list(x = x, I, value = value))
		)
	    }
	    ## Map.
	    if (is.double(I))
		I <- trunc(I)
	    if (any(I < 0L))
		stop("negative values are not allowed in a matrix subscript")
	    ## Remove rows with zero subscripts.
	    I <- I[.Call(R_all_row, I > 0L, FALSE),, drop = FALSE]
	    if (!nrow(I))
		return(x)
	    ## NOTE NAs cannot be introduced by coercion as
	    ##      long as the bounds are integer.
	    if (any(I > rep(x$dim, each = nrow(I))))
		stop("subscript out of bounds")
	    storage.mode(I) <- "integer"
    } else {
	if (na != nd + 2L)
	    stop("incorrect number of dimensions")
        ## Replace missing dimensions. 
        args <- substitute(list(...))
	for (k in seq_along(args)[-1L])
	    if (identical(as.character(args[[k]]), ""))
		if (.protect)
		    stop("Missing dimensions disabled for this object.")
		else
		    args[[k]] <- seq_len(x$dim[k - 1L])
        args <- eval(args)
	if (!all(sapply(args, is.numeric)))
	    stop("Only numeric subscripting is implemented.")
	## Replace negative subscripts.
	for (k in seq_along(args)) {
	    ## Map.
	    if (is.double(args[[k]]))
		args[[k]] <- trunc(args[[k]]) 
	    if (.protect) {
		if (any(args[[k]] < 0L))
		    stop("Negative subscripting disabled for this object.")
	    } else
		if (all(args[[k]] <= 0L))
		    args[[k]] <- seq_len(x$dim[k])[args[[k]]]
		else
		    if (!all(args[[k]] >= 0L))
			stop("only 0's may be mixed with negative subscripts")
	}
	## Expand.
	args <- matrix(
	    unlist(expand.grid(args), use.names = FALSE),
	    ncol = length(args)
	)
	return(
	    do.call("[<-.simple_sparse_array", 
		    list(x = x, args, value = value))
	)
    }

    ## Recycling.
    if (nrow(I) %% length(value))
	warning("number of items to replace is not a multiple of replacement length")
    V <- rep(value, length.out = nrow(I))

    ## Merge. 
    ##
    ## Emulates subsequent assignments of a sequence
    ## of replacement values with duplicate cell 
    ## indexes.
    I <- rbind(x$i, I)
    k <- .Call(R_match_matrix, I, NULL, NULL)
    k <- !duplicated(k[[1L]], fromLast = TRUE)
    I <- I[k,, drop = FALSE]
    V <- c(x$v, V)[k]

    ## Remove ZERO entries.
    k <- V == vector(typeof(V), 1L)
    if (any(k)) {
	k <- !k
	I <- I[k,, drop = FALSE]
	V <- V[k]
    }

    simple_sparse_array(
	v = V,
	i = I,
	dim = x$dim,
	dimnames = x$dimnames
    )
}

##
`[<-.simple_triplet_matrix` <- 
function(x, i, j, value) {
    x <- structure(list(
	    v = x$v, 
	    i = cbind(x$i, x$j), 
	    dim = c(x$nrow, x$ncol),
	    dimnames = x$dimnames
	),
	class =  "simple_sparse_array"
    )
    x <- list(x = x)
    if (!missing(i))
	x$i <- i
    if (!missing(j))
	x$j <- j
    x$value <- value
    x <- do.call("[<-.simple_sparse_array", x)
    if (inherits(x, "simple_sparse_array")) {
	x$i <- .Call(R_split_col, x$i)
	x <- structure(list(
		v = x$v, 
		i = x$i[[1L]], 
		j = x$i[[2L]], 
		nrow = x$dim[1L], 
		ncol = x$dim[2L], 
		dimnames = x$dimnames
	    ),
	    class = "simple_triplet_matrix"
	)
    }
    x
}

###
