## A simple class for sparse (triplet) matrices.

## Mostly intended for being able to take advantage of LP solvers which
## allow for sparse specifictions of (possible rather large) constraint
## matrices.

simple_triplet_matrix <-
function(i, j, v, nrow = max(i), ncol = max(j), dimnames = NULL)
{
    structure(list(i = as.integer(i), j = as.integer(j), v = v,
                   nrow = as.integer(nrow), ncol = as.integer(ncol),
                   dimnames = dimnames),
              class = "simple_triplet_matrix")
}

as.simple_triplet_matrix <-
function(x)
    UseMethod("as.simple_triplet_matrix")

as.simple_triplet_matrix.simple_triplet_matrix <- identity

as.simple_triplet_matrix.matrix <-
function(x)
{
    if(!prod(dim(x)))
        return(simple_triplet_matrix(integer(), integer(), c(x),
                                     nrow = nrow(x), ncol = ncol(x),
                                     dimnames = dimnames(x)))
    ind <- which(is.na(x) | (x != vector(typeof(x), 1L)), arr.ind = TRUE)
    dimnames(ind) <- NULL
    simple_triplet_matrix(ind[, 1L], ind[, 2L], x[ind],
                          nrow = nrow(x), ncol = ncol(x),
                          dimnames = dimnames(x))
}

as.simple_triplet_matrix.dgTMatrix <-
function(x)
{
    simple_triplet_matrix(x@i + 1L, x@j + 1L, x@x,
                          x@Dim[1L], x@Dim[2L], x@Dimnames)
}

## We could also simply write a method to coerce to a dgTMatrix, based
## on something like
##  new("dgTMatrix",
##       i = as.integer(i - 1),
##       j = as.integer(j - 1),
##       x = v,
##       Dim = c(nrow, ncol))
## (Note that these have C-style indices starting at zero.)

as.matrix.simple_triplet_matrix <-
function(x, ...)
{
    nr <- x$nrow
    nc <- x$ncol
    y <- matrix(vector(typeof(x$v), nr * nc), nr, nc)
    y[cbind(x$i, x$j)] <- x$v
    dimnames(y) <- x$dimnames
    y
}

is.simple_triplet_matrix <-
function(x)
    inherits(x, "simple_triplet_matrix")

is.numeric.simple_triplet_matrix <-
function(x)
    is.numeric(x$v)


Ops.simple_triplet_matrix <-
function(e1, e2)
{
    ## Currently, we only implement the following (for numeric
    ## operands):

    ## * Unary plus and minus.
    ##
    ## * Addition, subtraction and multiplication of two compatible
    ##   simple triplet matrices (or operands coercible to these).
    ##   [Division by a simple triplet matrix typically involves
    ##   division by zero and hence is not provided.]
    ##
    ## * Multiplication and division of a simple triplet matrix x by a
    ##   number or a vector of length nrow(x) (allowing to conveniently
    ##   scale the rows of a numeric simple triplet matrix).
    ##
    ## * Non-equality comparison of a simple triplet matrix with 0.
    ##
    ## * Comparisons of the elements of a simple triplet matrix with a
    ##   number.
    ##
    ## More could be added (but note that the elements could have
    ## arbitrary modes).

    op <- as.character(.Generic)

    if(nargs() == 1L) {
        if(op == "+") return(e1)
        if(op == "-") {
            e1$v <- - e1$v
            return(e1)
        }
        stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    }

    if(!(op %in% c("+", "-", "*", "/", "^",
                   "==", "!=", "<", "<=", ">", ">=")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    ## Require numeric operands for the arithmetic operators.
    if(!is.numeric(e1) || !is.numeric(e2))
        stop("Not implemented.")

    if(op %in% c("==", "!=", "<", "<=", ">", ">=")) {
        if(length(v <- e2) == 1L) {
            if(is.na(v))
                stop("NA/NaN handling not implemented.")
            ind <- if(do.call(.Generic, list(0, v))) {
                ## This inverts the sparse storage advantage, and hence
                ## will typically be inefficient.  Need to find the row
                ## and column positions of the zero entries.
                m <- matrix(1, e1$nrow, e1$ncol)
                m[cbind(e1$i, e1$j)] <- 0
                which(m != 0, arr.ind = TRUE)
            } else integer()
            e1$v <- do.call(.Generic, list(e1$v, v))
            pos <- which(!e1$v)
            if(length(pos)) {
                e1$i <- e1$i[-pos]
                e1$j <- e1$j[-pos]
                e1$v <- e1$v[-pos]
            }
            if(n <- NROW(ind)) {
                e1$i <- c(e1$i, ind[, 1L])
                e1$j <- c(e1$j, ind[, 2L])
                e1$v <- c(e1$v, rep.int(TRUE, n))
            }
            return(e1)
        }
        stop("Not implemented.")
    }

    if(op == "^") {
        ## Allow for taking (single) positive exponents.
        if(is.object(e2) || (length(e2) != 1L) ||
           !is.finite(e2) || (e2 <= 0))
            stop("Not implemented.")
        e1$v <- e1$v ^ e2
        return(e1)
    }

    .make_dimnames <- function(e1, e2) {
        if(is.null(rnms <- rownames(e1)))
            rnms <- rownames(e2)
        if(is.null(cnms <- colnames(e1)))
            cnms <- colnames(e2)
        if(is.null(rnms) && is.null(cnms))
            NULL
        else
            list(rnms, cnms)
    }

    .reduce <- function(x) {
	ind <- which(x$v == vector(typeof(x$v), 1L))
	if(length(ind)) {
	    x$i <- x$i[-ind]
            x$j <- x$j[-ind]
            x$v <- x$v[-ind]
        }
        x
    }

    ## Obviously, the following could be generalized ...

    if(op == "*") {
        if(!is.object(e1)) {
	    if(!all(is.finite(e1)))
		return(as.simple_triplet_matrix(e1 * as.matrix(e2)))
            if(length(e1) == 1L) {
                e2$v <- e2$v * e1
                return(.reduce(e2))
            }
            if(length(e1) == e2$nrow) {
                e2$v <- e2$v * e1[e2$i]
                return(.reduce(e2))
            }
	    if (is.matrix(e1)) {
		if (!all(dim(e1) == c(e2$nrow, e2$ncol)))
		    stop("Incompatible dimensions.")
		e2$v <- e2$v * e1[cbind(e2$i, e2$j)]
		return(.reduce(e2))
	    }
	    stop("Not implemented.")
        }
        if(!is.object(e2)) {
            if(!all(is.finite(e2)))
                return(as.simple_triplet_matrix(as.matrix(e1) * e2))
            if(length(e2) == 1L) {
                e1$v <- e1$v * e2
                return(.reduce(e1))
            }
            if(length(e2) == e1$nrow) {
                e1$v <- e1$v * e2[e1$i]
                return(.reduce(e1))
            }
	    if (is.matrix(e2)) {
		if (!all(dim(e2) == c(e1$nrow, e1$ncol)))
		    stop("Incompatible dimensions.")
		e1$v <- e1$v * e2[cbind(e1$i, e1$j)]
		return(.reduce(e1))
	    }
	    stop("Not implemented.")
        }
        ## This leaves multiplying two simple triplet matrices.
        e1 <- as.simple_triplet_matrix(e1)
        e2 <- as.simple_triplet_matrix(e2)
        ## Check dimensions: currently, no recycling.
        if(((nr <- e1$nrow) != e2$nrow) || ((nc <- e1$ncol) != e2$ncol))
            stop("Incompatible dimensions.")
        if(length(e1$v) < length(e2$v)) {
            ## Swap e1 and e2 so that duplicated indices can be found
            ## more efficiently.
            e3 <- e1
            e1 <- e2
            e2 <- e3
        }
        ## Find duplicated indices.
        pos <- match(paste(e2$i, e2$j, sep = "\r"),
                     paste(e1$i, e1$j, sep = "\r"),
                     nomatch = 0L)
        ind <- which(pos > 0L)
        if(!all(is.finite(e1$v)) || !all(is.finite(e2$v))) {
	    ## Augment and reduce
	    e2$i <- c(e2$i[ind], e2$i[-ind], e1$i[-pos])
	    e2$j <- c(e2$j[ind], e2$j[-ind], e1$j[-pos])
	    e2$v <- c(e2$v[ind] * e1$v[pos],
		vector(typeof(e2$v), 1L) * c(e2$v[-ind], e1$v[-pos]))
	    e2$dimnames <- .make_dimnames(e1, e2)
	    return(.reduce(e2))
	} else
	    return(simple_triplet_matrix(e2$i[ind], e2$j[ind],
                                         e2$v[ind] * e1$v[pos],
                                         nr, nc, .make_dimnames(e1, e2)))
    }

    ## This is slightly inefficent but special value handling is already
    ## in place.  Note v / 0 = v * 0^(-1) = v * Inf.
    if(op == "/") {
        if(!is.object(e2))
	    return(e1 * e2^(-1))
	e2 <- as.matrix(e2)
	if (!is.object(e1))
	    return(as.simple_triplet_matrix(e1 * e2^(-1)))
	return(e1 * e2^(-1))
    }

    ## This leaves adding and subtracting two simple triplet matrices.
    e1 <- as.simple_triplet_matrix(e1)
    e2 <- if(op == "+")
        as.simple_triplet_matrix(e2)
    else
        as.simple_triplet_matrix(-e2)
    ## Check dimensions: currently, no recycling.
    if((e1$nrow != e2$nrow) || (e1$ncol != e2$ncol))
        stop("Incompatible dimensions.")
    if(length(e1$v) < length(e2$v)) {
        ## Swap e1 and e2 so that duplicated indices can be found more
        ## efficiently.
        e3 <- e1
        e1 <- e2
        e2 <- e3
    }
    ## Find duplicated indices.
    pos <- match(paste(e2$i, e2$j, sep = "\r"),
                 paste(e1$i, e1$j, sep = "\r"),
                 nomatch = 0L)
    ind <- which(pos == 0L)
    ## Notice 0 + special value = special value.
    e1$v[pos] <- e1$v[pos] + e2$v[pos > 0L]
    e1$i <- c(e1$i, e2$i[ind])
    e1$j <- c(e1$j, e2$j[ind])
    e1$v <- c(e1$v, e2$v[ind])
    e1$dimnames <- .make_dimnames(e1, e2)
    e1
}

dim.simple_triplet_matrix <-
function(x)
    c(x$nrow, x$ncol)

## <TODO>
## Add a dim setter.
## </TODO>

dimnames.simple_triplet_matrix <-
function(x)
    x$dimnames

`dimnames<-.simple_triplet_matrix` <-
function(x, value)
{
    if(!is.null(value)) {
        ## Should be a list of length 2.
        if(!is.list(value) || (length(value) != 2L))
            stop("Invalid dimnames.")
        ind <- sapply(value, length) == 0L
        if(all(ind))
            value <- NULL
        else {
            dnx <- vector("list", 2L)
            dnx[!ind] <- lapply(value[!ind], as.character)
        }
    }
    if(is.null(value))
        x["dimnames"] <- list(NULL)
    else
        x$dimnames <- dnx
    x
}

`[.simple_triplet_matrix` <-
function(x, i, j, drop = FALSE)
{
    ## (Well, we certainly don't drop ...)

    ## (See e.g. `[.data.frame` for the trickeries of subscript methods:
    ## e.g., 
    ##   x[i = sample.int(nr, k), , drop = FALSE]
    ## counts 4 arguments (x, i, j and drop) where j is missing ...

    na <- nargs() - !missing(drop)
    if ((na == 1L) || (na == 2L) && missing(i)) 
        return(x)
    
    nr <- x$nrow
    nc <- x$ncol

    if(na == 2L) {
        ## Single index subscripting.
        ## Mimic subscripting matrices: no named argument handling in
        ## this case.
        if(is.character(i))
            stop("Character subscripting currently not implemented.")
        if(!is.matrix(i)) {
            if(is.logical(i))
                i <- which(rep(i, length.out = nr))
            ## Let's hope we have a vector.
            ## What if we have both negatives and positives?
            if(all(i >= 0)) {
                i <- i[i > 0]
                out <- vector(mode = typeof(x$v), length = length(i))
                pos <- match(i, (x$j - 1L) * nr + x$i, 0L)
                out[pos > 0L] <- x$v[pos]
            } else if(all(i <= 0)) {
                out <- vector(mode = typeof(x$v), nr * nc)
                out[(x$j - 1L) * nr + x$i] <- x$v
                out <- out[i]
            }
            else stop("Cannot mix positive and negative subscripts.")
        }
        else {
            ## Note that negative values are not allowed in a matrix
            ## subscript.
            if((ncol(i) != 2L) || (any(i < 0)))
                stop("Invalid subscript.")
            ## Rows containing zero indices can be dropped.
            ## Rows with NA indices should give NA (at least for
            ## non-recursive x).
            i <- i[!apply(i == 0, 1L, any), , drop = FALSE]
            out <- vector(mode = typeof(x$v), length = nrow(i))
            ##  pi <- match(i[, 1L], x$i)
            ##  pj <- match(i[, 2L], x$j)
            ##  ind <- which(pi == pj)
            ##  out[ind] <- x$v[pi[ind]]
            pos <- match(paste(i[,1L], i[,2L], sep = "\r"),
                         paste(x$i, x$j, sep = "\r"),
                         nomatch = 0L)
            out[pos > 0L] <- x$v[pos]
        }
    }
    else {
        ## Two index subscripting is rather tricky, as it can also be
        ## used for rearranging and "recycling" rows and columns.  Let
        ## us not support the latter for now, so that selected rows and
        ## columns must be unique.
        if(missing(i)) {
            pos <- rep.int(TRUE, length(x$v))
            pi <- seq_len(nr)
        }
        else {
            if(is.logical(i))
                i <- which(rep(i, length.out = nr))
            if(!is.numeric(i))
                stop("Two-index subscripting needs numeric or logical subscripts.")                
            pi <- seq_len(nr)
            if(all(i >= 0)) {
                i <- i[i > 0]
                if(any(duplicated(i)))
                    stop("Repeated indices currently not allowed.")
            } else if(all(i <= 0))
                i <- pi[i]
            else
                stop("Cannot mix positive and negative subscripts.")
            nr <- length(i)
            pos <- match(x$i, i, 0L) > 0L
            pi[i] <- seq_len(nr)
        }
        if(missing(j)) {
            pj <- seq_len(nc)
        }
        else {
            if(is.logical(j))
                j <- which(rep(j, length.out = nc))
            if(!is.numeric(j))
                stop("Two-index subscripting needs numeric or logical subscripts.")                
            pj <- seq_len(nc)
            if(all(j >= 0)) {
                j <- j[j > 0]
                if(any(duplicated(j)))
                    stop("Repeated indices currently not allowed.")
            } else if(all(j <= 0))
                j <- pj[j]
            else
                stop("Cannot mix positive and negative subscripts.")
            nc <- length(j)
            pos <- pos & (match(x$j, j, 0L) > 0L)
            pj[j] <- seq_len(nc)
        }
        if(!is.null(dnx <- x$dimnames))
            dnx <- list(dnx[[1L]][i], dnx[[2L]][j])

        out <- simple_triplet_matrix(pi[x$i[pos]], pj[x$j[pos]],
                                     x$v[pos], nr, nc, dnx)
    }

    out
}

rbind.simple_triplet_matrix <-
function(..., deparse.level = 1L)
{
    args <- lapply(Filter(Negate(is.null), list(...)),
                   as.simple_triplet_matrix)
    ## Ignore 'deparse.level' ...
    out <- Reduce(function(x, y) {
        if((nc <- ncol(x)) != ncol(y))
            stop("Numbers of columns of matrices must match.")
        nr <- nrow(x)
        simple_triplet_matrix(c(x$i, y$i + nr),
                              c(x$j, y$j),
                              c(x$v, y$v),
                              nrow = nr + nrow(y), ncol = nc)
    }, args)
    ## Handle dimnames in one final step.
    rnms <- lapply(args, rownames)
    rnms <- if(!any(sapply(rnms, is.null)))
        do.call("c", rnms)
    else
        NULL
    cnms <- Find(Negate(is.null), lapply(args, colnames))
    dimnames(out) <- list(rnms, cnms)
    out
}

cbind.simple_triplet_matrix <-
function(..., deparse.level = 1L)
{
    args <- lapply(Filter(Negate(is.null), list(...)),
                   as.simple_triplet_matrix)
    ## Ignore 'deparse.level' ...
    out <- Reduce(function(x, y) {
        if((nr <- nrow(x)) != nrow(y))
            stop("Numbers of rows of matrices must match.")
        nc <- ncol(x)
        simple_triplet_matrix(c(x$i, y$i),
                              c(x$j, y$j + nc),
                              c(x$v, y$v),
                              nrow = nr, ncol = nc + ncol(y))
    }, args)
    ## Handle dimnames in one final step.
    rnms <- Find(Negate(is.null), lapply(args, rownames))
    cnms <- lapply(args, rownames)
    cnms <- if(!any(sapply(cnms, is.null)))
        do.call("c", cnms)
    else
        NULL
    dimnames(out) <- list(rnms, cnms)
    out
}

t.simple_triplet_matrix <-
function(x)
    simple_triplet_matrix(x$j, x$i, x$v, x$ncol, x$nrow, rev(x$dimnames))

duplicated.simple_triplet_matrix <-
function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
    ## We could use the duplicated method for class matrix, but at the
    ## expense of going from sparse to dense ...
    if(!is.logical(incomparables) || incomparables)
        .NotYetUsed("incomparables != FALSE")
    if(MARGIN == 1L) {
        i <- x$i
        j <- x$j
        len <- x$nrow
    } else if(MARGIN == 2L) {
        i <- x$j
        j <- x$i
        len <- x$ncol
    } else
        stop("Invalid margin.")
    o <- order(i, j)
    y <- split(paste(j[o], x$v[o], sep = "\r"), i[o])
    tmp <- character(len)
    names(tmp) <- seq_along(tmp)
    tmp[names(y)] <- sapply(y, paste, collapse = "\r")
    duplicated(tmp, fromLast = fromLast)
}

unique.simple_triplet_matrix <-
function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
    if(!is.logical(incomparables) || incomparables)
        .NotYetUsed("incomparables != FALSE")
    ind <- !duplicated(x, MARGIN = MARGIN, fromLast = fromLast)
    if(MARGIN == 1L)
        x[which(ind), ]
    else
        x[, which(ind)]
}

c.simple_triplet_matrix <-
function(..., recursive = FALSE)
{
    args <- list(...)
    ind <- which(sapply(args, inherits, "simple_triplet_matrix"))
    args[ind] <-
        lapply(args[ind],
               function(x) {
                   y <- vector(typeof(x$v), x$nrow * x$ncol)
                   y[x$i + (x$j - 1L) * x$nrow] <- x$v
                   y
               })
    do.call("c", args)
}

print.simple_triplet_matrix <-
function(x, ...)
{
    writeLines(sprintf("A %s simple triplet matrix.",
                       paste(dim(x), collapse = "x")))
    invisible(x)
}

## Utitilies for creating special simple triplet matrices:

simple_triplet_zero_matrix <-
function(nrow, ncol = nrow, mode = "double")
    simple_triplet_matrix(integer(), integer(), vector(mode, 0L),
                          nrow, ncol)

simple_triplet_diag_matrix <-
function(v, nrow = length(v))
{
    v <- rep(v, length.out = nrow)
    i <- seq_len(nrow)
    simple_triplet_matrix(i, i, v, nrow, nrow)
}

