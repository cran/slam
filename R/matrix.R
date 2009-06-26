## A simple class for sparse (triplet) matrices.

## Mostly intended for being able to take advantage of LP solvers which
## allow for sparse specifictions of (possible rather large) constraint
## matrices.

simple_triplet_matrix <-
function(i, j, v, nrow = max(i), ncol = max(j), dimnames = NULL)
{
    structure(list(i = i, j = j, v = v, nrow = nrow, ncol = ncol,
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

## We could also simply write a method to coerce to a dgTMatrix, based
## on something like
##  new("dgTMatrix",
##       i = as.integer(i - 1),
##       j = as.integer(j - 1),
##       x = v,
##       Dim = c(nrow, ncol))
## (Note that these have C-style indices starting at zero.)

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

    if(!(op %in% c("+", "-", "*", "/",
                   "==", "!=", "<", "<=", ">", ">=")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    ## Require numeric operands for the arithmetic operators.
    if(!is.numeric(e1) || !is.numeric(e2))
        stop("Not implemented.")

    if(op %in% c("==", "!=", "<", "<=", ">", ">=")) {
        if(length(v <- e2) == 1L) {
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

    ## Obviously, the following could be generalized ...

    if(op == "*") {
        if(!is.object(e1)) {
            if(length(e1) == 1L) {
                e2$v <- e2$v * e1
                return(e2)
            }
            if(length(e1) == e2$nrow) {
                e2$v <- e2$v * e1[e2$i]
                return(e2)
            }
            stop("Not implemented.")
        }
        if(!is.object(e2)) {
            if(length(e2) == 1L) {
                e1$v <- e1$v * e2
                return(e1)
            }
            if(length(e2) == e1$nrow) {
                e1$v <- e1$v * e2[e1$i]
                return(e1)
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
        return(simple_triplet_matrix(e2$i[ind], e2$j[ind],
                                     e2$v[ind] * e1$v[pos],
                                     nr, nc, .make_dimnames(e1, e2)))
    }

    if(op == "/") {
        if(!is.object(e2)) {
            if(length(e2) == 1L) {
                e1$v <- e1$v / e2
                return(e1)
            }
            if(length(e2) == e1$nrow) {
                e1$v <- e1$v / e2[e1$i]
                return(e1)
            }
        }
        stop("Not implemented.")
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
function(x, i, j, ...)
{
    ## (Well, we certainly don't drop ...)

    ## Note that calling x[] with a simple triplet matrix x will call
    ## the subscript method with args x and missing ...
    na <- nargs()
    if((na == 1L) || (na == 2L) && missing(i))
        return(x)

    nr <- x$nrow
    nc <- x$ncol

    if(na == 2L) {
        ## Single index subscripting.
        if(is.logical(i))
            stop("Logical subscripting currently not implemented.")
        else if(is.character(i))
            stop("Character subscripting currently not implemented.")
        else if(!is.matrix(i)) {
            ## Let's hope we have a vector.
            ## What if we have both negatives and positives?
            if(all(i >= 0)) {
                i <- i[i > 0]
                out <- vector(mode = typeof(x$v), length = length(i))
                pos <- match(i, (x$j - 1L) * nr + x$i, 0)
                out[pos > 0] <- x$v[pos]
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
            i <- i[!apply(i == 0, 1L, any), , drop = FALSE]
            out <- vector(mode = typeof(x$v), length = nrow(i))
            ##  pi <- match(i[, 1L], x$i)
            ##  pj <- match(i[, 2L], x$j)
            ##  ind <- which(pi == pj)
            ##  out[ind] <- x$v[pi[ind]]
            pos <- match(paste(i[,1L], i[,2L], sep = "\r"),
                         paste(x$i, x$j, sep = "\r"),
                         nomatch = 0L)
            out[pos > 0] <- x$v[pos]
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
        else if(!is.numeric(i))
            stop("Only numeric two-index subscripting is implemented.")
        else {
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
            pos <- match(x$i, i, 0) > 0
            pi[i] <- seq_len(nr)
        }
        if(missing(j)) {
            pj <- seq_len(nc)
        }
        else if(!is.numeric(j))
            stop("Only numeric two-index subscripting is implemented.")
        else {
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
            pos <- pos & (match(x$j, j, 0) > 0)
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
    do.call(c, args)
}

print.simple_triplet_matrix <-
function(x, ...)
{
    writeLines(sprintf("\nA %s simple triplet matrix.\n",
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

