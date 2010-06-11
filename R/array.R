## A simple class for sparse arrays.

## Not very useful yet: need at least a subscript method.
## (Unfortunately, additional methods such as for rowSums/colSums or
## apply, etc., are not straightforward to add in an S3 world ...)

simple_sparse_array <-
function(i, v, dim = NULL, dimnames = NULL)
{
    ## <FIXME>
    ## Add some sanity checking eventually ...
    ## i should be a matrix of indices (non-"zero" entries).
    ## v should be a "vector" of non-zero values, with length equal to
    ## the number of rows of i.
    ## </FIXME>
    if(is.null(dim)) dim <- apply(i, 2L, max)
    ## <FIXME>
    ## Add checks for dimnames: should be NULL or a list of entries
    ## which are either NULL or character vectors as long as the
    ## corresponding dim.
    ## </FIXME>
    ssa <- list(i = i, v = v, dim = dim, dimnames = dimnames)
    class(ssa) <- "simple_sparse_array"
    ssa
}

as.simple_sparse_array <-
function(x)
    UseMethod("as.simple_sparse_array")

as.simple_sparse_array.simple_sparse_array <- identity

as.simple_sparse_array.array <-
function(x)
{
    if(!prod(dim(x)))
        simple_sparse_array(array(integer(), dim(x)), c(x),
                            dim(x), dimnames(x))
    ind <- which(is.na(x) | (x != vector(typeof(x), 1L)), arr.ind = TRUE)
    dimnames(ind) <- NULL
    simple_sparse_array(ind, x[ind], dim(x), dimnames(x))
}

as.simple_sparse_array.default <-
function(x)
    as.simple_sparse_array(as.array(x))

as.array.simple_sparse_array <-
function(x, ...)
{
    v <- x$v
    dim <- x$dim
    y <- array(vector(typeof(v), prod(dim)), dim = dim,
               dimnames = x$dimnames)
    y[x$i] <- v
    y
}

is.simple_sparse_array <-
function(x)
    inherits(x, "simple_sparse_array")

Summary.simple_sparse_array <-
function(..., na.rm = FALSE)
{
    v <- unlist(lapply(list(...),
                       function(e) {
                           v <- as.simple_sparse_array(e)$v
                           if(length(v) < prod(dim(e)))
                               v <- c(v, vector(typeof(v), 1L))
                           v
                       }),
                recursive = FALSE)
    do.call(.Generic, list(v, na.rm = na.rm))
}

dim.simple_sparse_array <-
function(x)
    x$dim

dimnames.simple_sparse_array <-
function(x)
    x$dimnames

## <TODO>
## Add dim and dimnames setters.
## </TODO>

`[.simple_sparse_array` <-
function(x, ...)
{
    ## Note that calling x[] with a simple sparse array x will call the
    ## subscript method with args x and missing ...
    na <- nargs()
    if((na == 1L) || (na == 2L) && missing(..1))
        return(x)

    nd <- length(x$dim)
    spos <- function(i) {
        ## Scalar positions of array index matrices i in the usual row
        ## major ordering of arrays.
        cpd <- cumprod(x$dim)
        1L + row_sums((i - 1L) * rep(c(1L, cpd[-nd]), each = cpd[nd]))
    }

    if(na == 2L) {
        i <- ..1
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
                pos <- match(i, spos(x$i), 0L)
                out[pos > 0L] <- x$v[pos]
            } else if(all(i <= 0)) {
                out <- vector(mode = typeof(x$v), prod(x$dim))
                out[spos(x$i)] <- x$v
                out <- out[i]
            }
            else stop("Cannot mix positive and negative subscripts.")
        }
        else {
            ## Note that negative values are not allowed in a matrix
            ## subscript.
            if((ncol(i) != nd) || (any(i < 0)))
                stop("Invalid subscript.")
            i <- i[!apply(i == 0, 1L, any), , drop = FALSE]
            out <- vector(mode = typeof(x$v), length = nrow(i))
            ## This is not really the fastest way to match rows, but is
            ## there an obvious better one?
            pos <- match(split(i, row(i)), split(x$i, row(x$i)), 0L)
            out[pos > 0L] <- x$v[pos]
        }
    }
    else {
        if(na != (nd + 1L))
            stop("Incorrect number of dimensions.")
        ## Figure out the missing arguments (if any).
        args <- substitute(list(...))
        ## Replace missing arguments by NULL for now.
        args[sapply(args,
                    function(a)
                    (length(a) == 1L) && (as.character(a) == ""))] <-
                        list(NULL)
        ## (Could also test args for identical(as.character(a), "").)
        ## And evaluate.
        args <- eval(args)
        ## Ready to go.
        dx <- x$dim
        pos <- rep.int(TRUE, length(x$v))
        ind <- lapply(dx, seq_len)
        for(k in seq_len(nd)) {
            i <- args[[k]]              # Given indices.
            if(is.null(i)) next
            else if(!is.numeric(i))
                stop("Only numeric multi-index subscripting is implemented.")
            else {
                if(all(i >= 0)) {
                    i <- i[i > 0]
                    if(any(duplicated(i)))
                        stop("Repeated indices currently not allowed.")
                } else if(all(i <= 0))
                    i <- seq_len(dx[k])[i]
                else
                    stop("Cannot mix positive and negative subscripts.")
                ind[[k]] <- i
                dx[k] <- length(i)
                j <- match(x$i[, k], i, 0L)
                x$i[j > 0L, k] <- seq_along(i)[j]
                pos <- pos & (j > 0L)
            }
        }
        if(!is.null(dnx <- x$dimnames))
            dnx[] <- Map("[", dnx, ind)
        out <- simple_sparse_array(x$i[pos, , drop = FALSE], x$v[pos],
                                   dx, dnx)
    }

    out

}

## <TODO>
## Add duplicated and unique methods for simple sparse arrays along the
## lines of the corresponding methods for simple triplet matrices.
## </TODO>

print.simple_sparse_array <-
function(x, ...)
{
    writeLines(sprintf("A %s simple sparse array.",
                       paste(dim(x), collapse = "x")))
    invisible(x)
}

mean.simple_sparse_array <-
function(x, ...)
{
    sum(x$v) / prod(dim(x))
}
