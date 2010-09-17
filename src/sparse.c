#include <R.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <time.h>

// ceeboo 2009/5,10,12 2010/1,5,6 
//

// test validity of list components.
int _valid_stm(SEXP x) {
    SEXP s = getAttrib(x, R_NamesSymbol);
    int ok = 
	strcmp(CHAR(STRING_ELT(s, 0)), "i") ||
	strcmp(CHAR(STRING_ELT(s, 1)), "j") ||
	strcmp(CHAR(STRING_ELT(s, 2)), "v") ||
	strcmp(CHAR(STRING_ELT(s, 3)), "nrow") ||
	strcmp(CHAR(STRING_ELT(s, 4)), "ncol") ||
    (LENGTH(s) > 5) ?
	strcmp(CHAR(STRING_ELT(s, 5)), "dimnames") : 0;
    if (ok) {
	s = VECTOR_ELT(x, 0);
	if (LENGTH(s) != LENGTH(VECTOR_ELT(x, 1)) ||
	    LENGTH(s) != LENGTH(VECTOR_ELT(x, 2)))
	    error("'i, j, v' different lengths");
	if (LENGTH(VECTOR_ELT(x, 3)) != 1 ||
	    LENGTH(VECTOR_ELT(x, 4)) != 1)
	    error("'nrow, ncol' invalid length");
	int *xi, *xj, nr, nc;
	xi = INTEGER(s);
	xj = INTEGER(VECTOR_ELT(x, 1));
	nr = INTEGER(VECTOR_ELT(x, 3))[0];
	nc = INTEGER(VECTOR_ELT(x, 4))[0];
	for (int k = 0; k < LENGTH(s); k++)
	    if (xi[k] < 1 || xi[k] > nr ||
		xj[k] < 1 || xj[k] > nc)
		error("'i, j' invalid");
	if (LENGTH(x) > 5) {
	    s = VECTOR_ELT(x, 5);
	    if (!isNull(s)) {
		if (LENGTH(s) != 2)
		    error("'dimnames' invalid length");
		if ((!isNull(VECTOR_ELT(s, 0)) &&
		      LENGTH(VECTOR_ELT(s, 0)) != nr) ||	
		    (!isNull(VECTOR_ELT(s, 1)) &&
		      LENGTH(VECTOR_ELT(s, 1)) != nc))
		    error("rownames, colnames invalid length'");
	    }
	}
    }
    return ok;
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
//          bailout function with y possibly coerced to REALSXP.
//       4) pkgEnv = NULL deactivates the bailout.
//       5) transpose 
//
SEXP tcrossprod_stm_matrix(SEXP x, SEXP R_y, SEXP pkgEnv, SEXP R_verbose,
							  SEXP R_transpose) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class simple_triplet_matrix");
    SEXP y = R_y;
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
    // coercing is in general not storage efficient, and therefore
    // assumes that y is not too large. on the other hand, as the
    // entries of y could be accessed multiple times, casting would
    // not be runtime efficient. if memory footprint is of concern
    // then the program flow should be further switch(ed).
    if (TYPEOF(y) != REALSXP)
	y = PROTECT(coerceVector(y, REALSXP));

    // check for special values
    SEXP r;
    double *_y = REAL(y);
    for (double *k = _y + LENGTH(y); _y < k; _y++)
	if (!R_FINITE(*_y)) {
	    if (isNull(pkgEnv))
		error("NA/NaN handling deactivated");
bailout:
	    r = eval(PROTECT(LCONS(install(".tcrossprod.bailout"),
			     LCONS(x,
			     LCONS(y, 
			     LCONS((R_transpose && *LOGICAL(R_transpose)) ?
				    R_transpose : ScalarLogical(FALSE),
				    R_NilValue))))), pkgEnv);
	    UNPROTECT(1);
	    if (y != R_y)
		UNPROTECT(1);
	    return r;
	}
    _y = REAL(y) - m;

    r = PROTECT(allocVector(REALSXP, n * m));
    memset(REAL(r), 0, sizeof(double) * n * m);
    double *_r = REAL(r) - m;

    int *_i, *_j;

    _i = INTEGER(VECTOR_ELT(x, 0));
    _j = INTEGER(VECTOR_ELT(x, 1));

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
    if (!R_transpose || !*LOGICAL(R_transpose)) {
	 v = r;
	_y = REAL(v);
	r = PROTECT(allocMatrix(REALSXP, n, m));
	_r = REAL(r);
	for (int i = 0; i < n * m; i++)
	    _r[i] = _y[i / n + (i % n) * m];
	UNPROTECT_PTR(v);
    } else {
	// NOTE we rely on setAttrib to not check if the dimnames
	//      are consistent with dim. 
	SEXP d = PROTECT(allocVector(INTSXP, 2));
	INTEGER(d)[0] = m;
	INTEGER(d)[1] = n;
	setAttrib(r, R_DimSymbol, d);
	UNPROTECT(1);
    }
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
    // swap dimnames
    if (R_transpose && *LOGICAL(R_transpose)) {
	dn = getAttrib(r, R_DimNamesSymbol);
	if (!isNull(dn)) {
	    SEXP t;
	    t = VECTOR_ELT(dn, 0);
	    SET_VECTOR_ELT(dn, 0, VECTOR_ELT(dn, 1));
	    SET_VECTOR_ELT(dn, 1, t);
	    dn = getAttrib(dn, R_NamesSymbol);
	    if (!isNull(dn)) {
		t = STRING_ELT(dn, 0);
		SET_STRING_ELT(dn, 0, STRING_ELT(dn, 1));
		SET_STRING_ELT(dn, 1, t);
	    }
	}
    }
#ifdef _TIME_H
    t3 = clock();
    if (R_verbose && *LOGICAL(R_verbose))
	Rprintf("tcrossprod_stm_matrix: %.3fs [%.3fs/%.3fs/%.3fs]\n", 
		((double) t3 - t0) / CLOCKS_PER_SEC,
		((double) t1 - t0) / CLOCKS_PER_SEC,
		((double) t2 - t1) / CLOCKS_PER_SEC,
		((double) t3 - t2) / CLOCKS_PER_SEC);
#endif
    UNPROTECT(1);
    if (y != R_y)
	UNPROTECT(1);

    return r;
}

// tcrossprod for some triplet matrices. 
//
// NOTES 1) y is not implemented.
//       2) pkgEnv = NULL deactivates the bailout to dense 
//          computation.
//
SEXP tcrossprod_stm_stm(SEXP x, SEXP y, SEXP pkgEnv, SEXP R_verbose) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class simple_triplet_matrix");
    if (!isNull(y))
	error("'y' not implemented");
    int *_i, *_j, *_n, k, f, n, m;
    double *_v, *_r;
    SEXP r, s;
#ifdef _TIME_H
    clock_t t2, t1, t0 = clock();
#endif
    s = VECTOR_ELT(x, 2);
    if (TYPEOF(s) != REALSXP)
	s = PROTECT(coerceVector(s, REALSXP));
    _v = REAL(s);
    for (k = 0; k < LENGTH(s); k++)
	if (!R_FINITE(_v[k])) {
	    if (isNull(pkgEnv))
		error("NA/NaN handling deactivated");
	    r = eval(PROTECT(LCONS(install(".tcrossprod.bailout"),
			     LCONS(x,
			     LCONS(y, 
			     LCONS(ScalarLogical(FALSE), 
				   R_NilValue))))), pkgEnv);
	    UNPROTECT(1);
	    if (s != VECTOR_ELT(x, 2))
		UNPROTECT(1);
	    return r;
	}

    n = INTEGER(VECTOR_ELT(x, 3))[0];
    if (!n) {
	if (s != VECTOR_ELT(x, 2))
	    UNPROTECT(1);
	return allocMatrix(REALSXP, 0, 0);
    }
    m = INTEGER(VECTOR_ELT(x, 4))[0];
    r = PROTECT(allocMatrix(REALSXP, n, n));
    memset(REAL(r), 0, sizeof(double) * n * n);
    if (LENGTH(x) > 5) {
	SEXP s = VECTOR_ELT(x, 5);
	if (!isNull(s)) {
	    SEXP t = VECTOR_ELT(s, 0);
	    if (!isNull(t)) {
		SEXP d;
		setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
		SET_VECTOR_ELT(d, 0, t);
		SET_VECTOR_ELT(d, 1, t);
		s = getAttrib(s, R_NamesSymbol);
		if (!isNull(s)) {
		    t = STRING_ELT(s, 0);
		    setAttrib(d, R_NamesSymbol, (s = allocVector(STRSXP, 2)));
		    SET_STRING_ELT(s, 0, t);
		    SET_STRING_ELT(s, 1, t);
		}
	    }
	}
    }
    if (!m || !LENGTH(s)) {
	UNPROTECT(1);
	if (s != VECTOR_ELT(x, 2))
	    UNPROTECT(1);
	return r;
    }
    // Arrange the data in blocks of equal column
    // indexes. Note that the order of and within
    // the blocks is not relevant (see below).

    _j = INTEGER(VECTOR_ELT(x, 1));	    // column indexes
    _n = INTEGER(PROTECT(allocVector(INTSXP, m + 1)));
    memset(_n, 0, sizeof(int) * (m + 1));
    for (k = 0; k < LENGTH(s); k++)
	_n[_j[k]]++;
    for (k = 1; k < m + 1; k++)
	_n[k] += _n[k-1];
    {
	int *__i;
	double *__v;

	__i = INTEGER(VECTOR_ELT(x, 0));    // row indexs
	__v = _v;

	_i = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(s))));
	_v = REAL(PROTECT(allocVector(REALSXP, LENGTH(s))));

	_n -= 1;
	for (k = 0; k < LENGTH(s); k++) {
	    int *__n = _n + _j[k];
	    _i[*__n] = __i[k];
	    _v[*__n] = __v[k];
	    (*__n)++;
	}
	// reset
	_n += 1;
	for (k = m; k > 0; k--)
	    _n[k] = _n[k-1];
	_n[0] = 0;
    }
#ifdef _TIME_H
    t1 = clock();
#endif
    // Aggregate the outer products of the columns.
    _r = REAL(r) - n - 1;
     f = _n[0];
    for (k = 1; k < m + 1; k++) {
	int l = _n[k];
	for (int j = f; j < l; j++) {
	    double  z =      _v[j],
		  *_z = _r + _i[j] * n;
	    for (int i = f; i < j + 1; i++)
		_z[_i[i]] += _v[i] * z;
	}
	f = l;
    }
    // Aggregate the lower and upper half.
    _r = REAL(r);
    for (k = 1; k < n; k++) {
	f = k * n;
	// NOTE the off-diagonal array indexes are i * n + k, 
	//	and k * n + i for i = 0, 1, ..., k-1. For the
	//      former (k - 1) * n + k < k * n  <=>  k < n,
	//      and adding k to the right-hand sides does
	//      change that.
	for (int i = k; i < f; i += n, f++) {
	    _r[f] += _r[i];
	    _r[i]  = _r[f];
	}
    }
#ifdef _TIME_H
    t2 = clock();
    if (R_verbose && *LOGICAL(R_verbose))
	Rprintf("_crossprod_stm: %.3fs [%.3fs/%.3fs]\n",
		((double) t2 - t0) / CLOCKS_PER_SEC,
		((double) t1 - t0) / CLOCKS_PER_SEC,
		((double) t2 - t1) / CLOCKS_PER_SEC);
#endif
    UNPROTECT(4);
    if (s != VECTOR_ELT(x, 2))
	UNPROTECT(1);

    return r;
}



//

