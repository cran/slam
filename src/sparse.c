#include <R.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <time.h>

// ceeboo 2009/5,10,12 2010/1 
//

// test order of list components
int _valid_stm(SEXP x) {
    x = getAttrib(x, R_NamesSymbol);
    return 
	strcmp(CHAR(STRING_ELT(x, 0)), "i") ||
	strcmp(CHAR(STRING_ELT(x, 1)), "j") ||
	strcmp(CHAR(STRING_ELT(x, 2)), "v") ||
	strcmp(CHAR(STRING_ELT(x, 3)), "nrow") ||
	strcmp(CHAR(STRING_ELT(x, 4)), "ncol") ||
    (LENGTH(x) > 5) ?
	strcmp(CHAR(STRING_ELT(x, 5)), "dimnames") : 0;
}

// row or column sums of some triplet matrix
//
SEXP _sums_stm(SEXP x, SEXP R_dim, SEXP R_na_rm) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class 'simple_triplet_matrix'");
    if (TYPEOF(R_dim) != INTSXP)
	error("'dim' not of type integer");
    if (TYPEOF(R_na_rm) != LGLSXP)
	error("'na.rm' not of type logical");

    int n, *i = NULL;
    
    switch ((n = *INTEGER(R_dim))) {
	case 1:
	    i = INTEGER(VECTOR_ELT(x, 0));
	    break;
	case 2:
	    i = INTEGER(VECTOR_ELT(x, 1));
	    break;
	default:
	    error("'dim' invalid");
    }
    n = INTEGER(VECTOR_ELT(x, n + 2))[0];

    // for the type of the return argument see the behavior
    // of rowSums and colSums for matrix.
    SEXP r = PROTECT(allocVector(REALSXP, n));

    memset(REAL(r), 0, sizeof(double) * n);
    // offset one-based indexing
    double *__r__ = REAL(r) - 1;

    SEXP _x_ = VECTOR_ELT(x, 2);

    switch (TYPEOF(_x_)) {
	case LGLSXP:
	case INTSXP: {
	    int v, *k, *__x__ = INTEGER(_x_);
	    if (*LOGICAL(R_na_rm)) {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    if ((v = *__x__) == NA_INTEGER)
			continue;
		    else 
			__r__[*i] += (double) v;
	    } else {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    __r__[*i] +=
		    // map NA
			 ((v = *__x__) == NA_INTEGER) ? NA_REAL : v;
	    }
	    break;
	}
	case REALSXP: {
	    double v, *k, *__x__ = REAL(_x_);
	    if (*LOGICAL(R_na_rm)) {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    if (ISNAN((v = *__x__)))
			continue;
		    else
			__r__[*i] += v;
	    } else
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    __r__[*i] += *__x__;
	    break;
	}
	default:
	    error("type of 'x' not supported");
    }

    SEXP d = (LENGTH(x) > 5) ? VECTOR_ELT(x, 5) : R_NilValue;
    if (!isNull(d)) {
	n = *INTEGER(R_dim);
	setAttrib(r, R_NamesSymbol, VECTOR_ELT(d, n - 1));
    }

    UNPROTECT(1);

    return r;
}

// tcrossprod for some triplet matrix and matrix
//
// NOTES 1) tcrossprod does not implement na.rm, so neither do we.
//       2) triplet on triplet does not fit in here.
//       3) if y = NULL or contains special values we call some
//          bailout function.
SEXP tcrossprod_stm_matrix(SEXP x, SEXP y, SEXP pkgEnv, SEXP R_verbose) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class simple_triplet_matrix");
    if (isNull(y))
	goto bailout;
    if (!isMatrix(y))
	error("'y' not of class matrix");

    int n, m;

    n = INTEGER(VECTOR_ELT(x, 4))[0];

    if (n != INTEGER(getAttrib(y, R_DimSymbol))[1])
	error("the number of columns of 'x' and 'y' do not conform");

    n = INTEGER(VECTOR_ELT(x, 3))[0];
    m = INTEGER(getAttrib(y, R_DimSymbol))[0];

#ifdef _TIME_H
    // code section times
    clock_t t3, t2, t1, t0 = clock();
#endif
    SEXP r = PROTECT(allocVector(REALSXP, n * m));

    memset(REAL(r), 0, sizeof(double) * n * m);
    double *_r = REAL(r) - m;

    int *_i, *_j;

    _i = INTEGER(VECTOR_ELT(x, 0));
    _j = INTEGER(VECTOR_ELT(x, 1));

    // coercing is in general not storage efficient, and therefore
    // assumes that y is not too large. on the other hand, as the
    // entries of y could be accessed multiple times, casting would
    // not be runtime efficient. if memory footprint is of concern
    // then the program flow should be further switch(ed).
    if (TYPEOF(y) != REALSXP)
	y = coerceVector(y, REALSXP);

    // check for special values
    double *_y = REAL(y);
    for (double *k = _y + LENGTH(y); _y < k; _y++)
	if (!R_FINITE(*_y)) {
	    UNPROTECT(1);
bailout:
	    return eval(LCONS(install(".tcrossprod.bailout"),
			LCONS(x,
			LCONS(y, R_NilValue))), pkgEnv);
	}
    _y = REAL(y) - m;

    // Notes 1) timings with Blas are better than without.
    //       2) For reasons not yet fully understood using
    //          a transposed result matrix is more runtime
    //          efficient.
    SEXP v = VECTOR_ELT(x, 2);
#ifdef _TIME_H
    t1 = clock();
#endif
    switch (TYPEOF(v)) {
	case LGLSXP:
	case INTSXP: {
	    int    *k, *__x = INTEGER(v);
	    double *l, *__r, *__y;
	    for (k = __x + LENGTH(v); __x < k; __x++, _i++, _j++) {
		__r = _r + *_i * m;
		__y = _y + *_j * m;

		for (l = __y + m; __y < l; __y++, __r++)
		    *__r += *__x * *__y;
	    }
	    break;
	}
	case REALSXP: {
	    double *k, *__x = REAL(v);
#ifdef R_BLAS_H
	    int l = 1, *_l = &l, *_m = &m;
#else
	    double *l, *__r, *__y;
#endif
	    for (k = __x + LENGTH(v); __x < k; __x++, _i++, _j++) {
#ifdef R_BLAS_H
		F77_NAME(daxpy)(_m, __x, _y + *_j * m, _l,
			                 _r + *_i * m, _l);
#else
		__r = _r + *_i * m;
		__y = _y + *_j * m;

		for (l = __y + m; __y < l; __y++, __r++)
		    *__r += *__x * *__y;
#endif
	    }
	    break;
	}
	default:
	    error("type of 'x' not supported");
    }
#ifdef _TIME_H
    t2 = clock();
#endif
    // transpose
     v = r;
    _y = REAL(v);
     r = PROTECT(allocMatrix(REALSXP, n, m));
    _r = REAL(r);
    for (int i = 0; i < n * m; i++)
	_r[i] = _y[i / n + (i % n) * m];
    UNPROTECT_PTR(v);
#ifdef _TIME_H
    t3 = clock();
    if (R_verbose && *LOGICAL(R_verbose))
	Rprintf("tcrossprod_stm_matrix: %.3fs [%.3fs/%.3fs/%.3fs]\n", 
		((double) t3 - t0) / CLOCKS_PER_SEC,
		((double) t1 - t0) / CLOCKS_PER_SEC,
		((double) t2 - t1) / CLOCKS_PER_SEC,
		((double) t3 - t2) / CLOCKS_PER_SEC);
#endif
    // set dimnames and names of dimnames.
    SEXP dn = (LENGTH(x) > 5) ? VECTOR_ELT(x, 5) : R_NilValue;

    if (!isNull(dn)) {
	SEXP d, dnn;

	dnn = getAttrib(dn, R_NamesSymbol);

	setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	SET_VECTOR_ELT(d, 0, VECTOR_ELT(dn, 0));

	dn = getAttrib(y, R_DimNamesSymbol);
	if (!isNull(dn)) {
	    SET_VECTOR_ELT(d, 1, VECTOR_ELT(dn, 0));

	    if (!isNull(dnn)) {
		SEXP t;
		setAttrib(d, R_NamesSymbol, (t = allocVector(STRSXP, 2)));
		SET_STRING_ELT(t, 0, STRING_ELT(dnn, 0));

		dnn = getAttrib(dn, R_NamesSymbol);
		if (!isNull(dnn))
		    SET_STRING_ELT(t, 1, STRING_ELT(dnn, 0));
		else
		    SET_STRING_ELT(t, 1, R_BlankString);
	    } else {
		dnn = getAttrib(dn, R_NamesSymbol);
		if (!isNull(dnn)) {
		    SEXP t;
		    setAttrib(d, R_NamesSymbol, (t = allocVector(STRSXP, 2)));
		    SET_STRING_ELT(t, 0, R_BlankString);
		    SET_STRING_ELT(t, 1, STRING_ELT(dnn, 0));
		}
	    }
	} else {
	    SET_VECTOR_ELT(d, 1, R_NilValue);

	    if (!isNull(dnn)) {
		SEXP t;
		setAttrib(d, R_NamesSymbol, (t = allocVector(STRSXP, 2)));
		SET_STRING_ELT(t, 0, STRING_ELT(dnn, 0));
		SET_STRING_ELT(t, 1, R_BlankString);
	    }
	}
    } else {
	dn = getAttrib(y, R_DimNamesSymbol);
	if (!isNull(dn)) {
	    SEXP d;
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, R_NilValue);
	    SET_VECTOR_ELT(d, 1, VECTOR_ELT(dn, 0));

	    dn = getAttrib(dn, R_NamesSymbol);
	    if (!isNull(dn)) {
		SEXP t;
		setAttrib(d, R_NamesSymbol, (t = allocVector(STRSXP, 2)));
		SET_STRING_ELT(t, 0, R_BlankString);
		SET_STRING_ELT(t, 1, STRING_ELT(dn, 0));
	    }
	}
    }
    UNPROTECT(1);

    return r;
}

//

