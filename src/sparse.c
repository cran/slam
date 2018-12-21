#include <R.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Complex.h>
#include <time.h>

// ceeboo 2009/5,10,12 2010/1,5,6 2011/2 2012/4,5 2013/10 2016/6
//

// remove attributes from payload vector (see src/main/coerce.c)
SEXP _unattr(SEXP x) {
    if (!isVector(x) || ATTRIB(x) == R_NilValue)
	return x;
    if (MAYBE_SHARED(x)) {
	SEXP s = x;
	SEXP a = PROTECT(ATTRIB(x));
	SET_ATTRIB(x, R_NilValue);
	x = duplicate(x);
	SET_ATTRIB(s, a);
	UNPROTECT(1);		/* a */
    } else
	SET_ATTRIB(x, R_NilValue);
    if (OBJECT(x))
	SET_OBJECT(x, 0);
    if (IS_S4_OBJECT(x))
	UNSET_S4_OBJECT(x);
    return x;
}

// test validity of payload vector
int _valid_v(SEXP x) {
    if (!isVector(x))
	error("'x' not a vector");
    int i;
    i = LENGTH(x);
    switch(TYPEOF(x)) {
	case LGLSXP:
	    // test for FALSE (see below)
	case INTSXP: 
	    {
		int *v = INTEGER(x);
		while (i-- > 0)
		    if (v[i] == 0)
			break;
	    }
	    break;
	case REALSXP: 
	    {
		double *v = REAL(x);
		while (i-- > 0)
		    if (v[i] == (double) 0)
			break;
	    }
	    break;
	case RAWSXP:
	    { 
		unsigned char *v = RAW(x);
		while (i-- > 0)
		    if (v[i] == (unsigned char) 0)
			break;
	    }
	    break;
	case CPLXSXP:
	    {
		Rcomplex *v = COMPLEX(x);
		while (i-- > 0)
		    if (v[i].i == (double) 0 && 
			v[i].r == (double) 0)
			break;
	    }
	    break;
	case EXPRSXP:
	case VECSXP:
	    while (i-- > 0)
		if (VECTOR_ELT(x, i) == R_NilValue)
		    break;
	    break;
	case STRSXP:
	    while (i-- > 0)
		if (STRING_ELT(x, i) == R_BlankString)
		    break;
	    break;
	default:
	    error("type of 'x' not implemented");
	
    }
    return i + 1;
}

// wrapper
SEXP __valid_v(SEXP x) {
    return ScalarLogical(_valid_v(x) == FALSE);
}


// test validity of list components.
int _valid_stm(SEXP x) {
    if (LENGTH(x) < 5)
	error("invalid number of components");
    SEXP s = getAttrib(x, R_NamesSymbol);
    int ok = 
	strcmp(CHAR(STRING_ELT(s, 0)), "i") ||
	strcmp(CHAR(STRING_ELT(s, 1)), "j") ||
	strcmp(CHAR(STRING_ELT(s, 2)), "v") ||
	strcmp(CHAR(STRING_ELT(s, 3)), "nrow") ||
	strcmp(CHAR(STRING_ELT(s, 4)), "ncol") ||
    ((LENGTH(s) > 5) ?
	strcmp(CHAR(STRING_ELT(s, 5)), "dimnames") : 0);
    if (!ok) {
	if (TYPEOF(VECTOR_ELT(x, 0)) != INTSXP ||
	    TYPEOF(VECTOR_ELT(x, 1)) != INTSXP ||
	    TYPEOF(VECTOR_ELT(x, 3)) != INTSXP ||
	    TYPEOF(VECTOR_ELT(x, 4)) != INTSXP)
	    error("'i, j, nrow, ncol' invalid type");
	if (!isVector(VECTOR_ELT(x, 2)))
	    error("'v' not a vector");
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
	if (nr < 0 || nr == NA_INTEGER ||
	    nc < 0 || nc == NA_INTEGER)
	    error("'nrow, ncol' invalid");
	for (int k = 0; k < LENGTH(s); k++)
	    if (xi[k] < 1 || xi[k] > nr ||
		xj[k] < 1 || xj[k] > nc)
		error("'i, j' invalid");
	if (LENGTH(x) > 5) {
	    s = VECTOR_ELT(x, 5);
	    if (!isNull(s)) {
		if (TYPEOF(s) != VECSXP)
		    error("'dimnames' invalid type");
		if (LENGTH(s) != 2)
		    error("'dimnames' invalid length");
		if ((!isNull(VECTOR_ELT(s, 0)) &&
		     (LENGTH(VECTOR_ELT(s, 0)) != nr ||
		   !isString(VECTOR_ELT(s, 0)))) ||	
		    (!isNull(VECTOR_ELT(s, 1)) &&
		     (LENGTH(VECTOR_ELT(s, 1)) != nc ||
		   !isString(VECTOR_ELT(s, 1)))))
		    error("'dimnames' component invalid length or type");
	    }
	}
    }
    return ok;
}

// wrapper
SEXP __valid_stm(SEXP x) {
    if (!inherits(x, "simple_triplet_matrix"))
	return ScalarLogical(FALSE);
    return ScalarLogical(_valid_stm(x) == FALSE);
}

// row or column sums of some triplet matrix
//
SEXP _sums_stm(SEXP x, SEXP R_dim, SEXP R_na_rm) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class 'simple_triplet_matrix'");
    if (TYPEOF(R_dim) != INTSXP)
	error("'dim' not of type integer");
    if (!LENGTH(R_dim))
	error("'dim' invalid length");
    if (TYPEOF(R_na_rm) != LGLSXP)
	error("'na.rm' not of type logical");
    if (!LENGTH(R_na_rm))
	error("'na.rm' invalid length");

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

    SEXP r = NULL;
    SEXP _x_ = VECTOR_ELT(x, 2);

    switch (TYPEOF(_x_)) {
	case LGLSXP:
	case INTSXP: {
	    // for the type of the return argument see the behavior
	    // of rowSums and colSums for matrix.
	    r = PROTECT(allocVector(REALSXP, n));

	    memset(REAL(r), 0, sizeof(double) * n);
	    // offset one-based indexing
	    double *__r__ = REAL(r) - 1;

	    int *k, *__x__ = INTEGER(_x_);
	    if (*LOGICAL(R_na_rm)) {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    if (*__x__ == NA_INTEGER)
			continue;
		    else 
			__r__[*i] += (double) *__x__;
	    } else {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    if (*__x__ == NA_INTEGER)
			__r__[*i]  = NA_REAL;	// map NA
		    else
			__r__[*i] += (double) *__x__;
	    }
	    break;
	}
	case REALSXP: {
	    r = PROTECT(allocVector(REALSXP, n));

	    memset(REAL(r), 0, sizeof(double) * n);
	    double *__r__ = REAL(r) - 1;

	    double *k, *__x__ = REAL(_x_);
	    if (*LOGICAL(R_na_rm)) {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    if (ISNAN(*__x__))
			continue;
		    else
			__r__[*i] += *__x__;
	    } else
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    __r__[*i] += *__x__;
	    break;
	}
	case CPLXSXP: {
	    r = PROTECT(allocVector(CPLXSXP, n));

	    memset(COMPLEX(r), 0, sizeof(Rcomplex) * n);
	    Rcomplex *__r__ = COMPLEX(r) - 1;

	    Rcomplex *k, *__x__ = COMPLEX(_x_);
	    if (*LOGICAL(R_na_rm)) {
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++)
		    if (ISNAN(__x__->r) || ISNAN(__x__->i))
			continue;
		    else {
			__r__[*i].r += __x__->r;
			__r__[*i].i += __x__->i;
		    }
	    } else
		for (k = __x__ + LENGTH(_x_); __x__ < k; __x__++, i++) {
		    __r__[*i].r += __x__->r;
		    __r__[*i].i += __x__->i;
		}

	    break;
	}
	default:
	    error("type of 'x' invalid");
    }

    SEXP d = (LENGTH(x) > 5) ? VECTOR_ELT(x, 5) : R_NilValue;
    if (!isNull(d)) {
	n = *INTEGER(R_dim);
	setAttrib(r, R_NamesSymbol, VECTOR_ELT(d, n - 1));
    }

    UNPROTECT(1);

    return r;
}

// tcrossprod for some triplet matrices. 
//
// NOTES 1) y is now implemented.
//       2) pkgEnv = NULL deactivates the bailout to dense 
//          computation.
//
SEXP tcrossprod_stm_stm(SEXP x, SEXP y, SEXP pkgEnv, SEXP R_verbose) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class simple_triplet_matrix");
    if (!isNull(y) &&
	(!inherits(y, "simple_triplet_matrix") || _valid_stm(y)))
	error("'y' not of class simple_triplet_matrix");
    int *_ix, *_jx, *_nx, k, fx, l, n, m;
    double *_vx, *_vy = NULL, *_r;
    SEXP r, vx, vy = NULL;

    l = INTEGER(VECTOR_ELT(x, 4))[0];
    if (!isNull(y) &&
	l != INTEGER(VECTOR_ELT(y, 4))[0])
	error("the number of columns of 'x' and 'y' do not conform");

#ifdef _TIME_H
    clock_t t2, t1, t0 = clock();
#endif
    vx = VECTOR_ELT(x, 2);
    if (TYPEOF(vx) != REALSXP)
	vx = PROTECT(coerceVector(vx, REALSXP));
    _vx = REAL(vx);
    for (k = 0; k < LENGTH(vx); k++)
	if (!R_FINITE(_vx[k])) {
	    if (isNull(pkgEnv))
		error("NA/NaN handling deactivated");
	    if (vx != VECTOR_ELT(x, 2))
		UNPROTECT(1);
	    r = eval(PROTECT(LCONS(install(".tcrossprod_bailout"),
		     PROTECT( CONS(x,
			      CONS(y, 
			      CONS(ScalarLogical(FALSE), 
				   R_NilValue)))))), pkgEnv);
	    UNPROTECT(2);
	    return r;
	}

    n = INTEGER(VECTOR_ELT(x, 3))[0];

    if (!isNull(y)) {
	vy = VECTOR_ELT(y, 2);
	if (TYPEOF(vy) != REALSXP)
	    vy = PROTECT(coerceVector(vy, REALSXP));
	_vy = REAL(vy);
	for (k = 0; k < LENGTH(vy); k++)
	    if (!R_FINITE(_vy[k])) {
		if (isNull(pkgEnv))
		    error("NA/NaN handling deactivated");
		if (vy != VECTOR_ELT(y, 2))
		    UNPROTECT(1);
		if (vx != VECTOR_ELT(x, 2))
		    UNPROTECT(1);
		r = eval(PROTECT(LCONS(install(".tcrossprod_bailout"),
			 PROTECT( CONS(x,
				  CONS(y, 
				  CONS(ScalarLogical(FALSE), 
				       R_NilValue)))))), pkgEnv);
		UNPROTECT(2);
		return r;
	    }

	m = INTEGER(VECTOR_ELT(y, 3))[0];
    } else
	m = n;

    r = PROTECT(allocMatrix(REALSXP, n, m));
    memset(REAL(r), 0, sizeof(double) * n * m);
    {
	SEXP sx, dx, sy, dy;
	sx = dx = sy = dy = R_NilValue; 
	if (LENGTH(x) > 5) {
	    sx = VECTOR_ELT(x, 5);
	    if (!isNull(sx)) {
		dx = VECTOR_ELT(sx, 0);
		sx = getAttrib(sx, R_NamesSymbol);
		if (!isNull(sx))
		    sx = STRING_ELT(sx, 0);
	    }
	}
	if (!isNull(y)) {
	    if (LENGTH(y) > 5) {
		sy = VECTOR_ELT(y, 5);
		if (!isNull(sy)) {
		    dy = VECTOR_ELT(sy, 0);
		    sy = getAttrib(sy, R_NamesSymbol);
		    if (!isNull(sy))
			sy = STRING_ELT(sy, 0);
		}
	    }
	} else {
	    sy = sx;
	    dy = dx;
	}
	if (!isNull(dx) || !isNull(dy)) {
	    SEXP d;
	    setAttrib(r, R_DimNamesSymbol, (d = allocVector(VECSXP, 2)));
	    SET_VECTOR_ELT(d, 0, dx);
	    SET_VECTOR_ELT(d, 1, dy);
	    if (!isNull(sx) || !isNull(sy)) {
		SEXP s;
		setAttrib(d, R_NamesSymbol, (s = allocVector(STRSXP, 2)));
		SET_STRING_ELT(s, 0, isNull(sx) ? R_BlankString : sx);
		SET_STRING_ELT(s, 1, isNull(sy) ? R_BlankString : sy);
	    }
	}
    }
    if (!l || !n || !LENGTH(vx) || 
	(!isNull(y) && (!m || !LENGTH(vy)))) {
	UNPROTECT(1);
	if (vx != VECTOR_ELT(x, 2))
	    UNPROTECT(1);
	if (!isNull(y) &&
	    vy != VECTOR_ELT(y, 2))
	    UNPROTECT(1);
	return r;
    }
    // Arrange the data in blocks of equal column
    // indexes. Note that the order within (of)
    // the blocks is not relevant (see below).

    _jx = INTEGER(VECTOR_ELT(x, 1));	    // column indexes
    _nx = INTEGER(PROTECT(allocVector(INTSXP, l + 1)));
    memset(_nx, 0, sizeof(int) * (l + 1));
    for (k = 0; k < LENGTH(vx); k++)
	_nx[_jx[k]]++;
    for (k = 1; k < l + 1; k++)
	_nx[k] += _nx[k-1];
    {
	int *__i;
	double *__v;

	__i = INTEGER(VECTOR_ELT(x, 0));    // row indexs
	__v = _vx;

	_ix = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(vx))));
	_vx = REAL(PROTECT(allocVector(REALSXP, LENGTH(vx))));

	_nx -= 1;
	for (k = 0; k < LENGTH(vx); k++) {
	    int *__n = _nx + _jx[k];
	    _ix[*__n] = __i[k];
	    _vx[*__n] = __v[k];
	    (*__n)++;
	}
	// reset
	_nx += 1;
	for (k = l; k > 0; k--)
	    _nx[k] = _nx[k-1];
	_nx[0] = 0;
    }

#ifdef _TIME_H
    t1 = clock();
#endif
    // Aggregate the outer products of the columns.
    if (isNull(y)) {
	_r = REAL(r) - n - 1;
	fx = _nx[0];
	for (k = 1; k < l + 1; k++) {
	    int lx = _nx[k];
	    for (int j = fx; j < lx; j++) {
		double  z =      _vx[j],
		      *_z = _r + _ix[j] * n;
		for (int i = fx; i < j + 1; i++)
		    _z[_ix[i]] += _vx[i] * z;
	    }
	    fx = lx;
	}
	// Aggregate the lower and upper half.
	_r = REAL(r);
	for (k = 1; k < n; k++) {
	    int j = k * n;
	    // NOTE the off-diagonal array indexes are i * n + k, 
	    //	    and k * n + i for i = 0, 1, ..., k-1. For the
	    //      former (k - 1) * n + k < k * n  <=>  k < n,
	    //      and adding k to the right sides does not
	    //      change that.
	    for (int i = k; i < j; i += n, j++) {
		_r[j] += _r[i];
		_r[i]  = _r[j];
	    }
	}
    } else {
	int *_iy, *_jy;

	_r = REAL(r) - n - 1;
	_iy = INTEGER(VECTOR_ELT(y, 0));
	_jy = INTEGER(VECTOR_ELT(y, 1));	    // column indexes
	for (k = 0; k < LENGTH(vy); k++) {
	    int    j = _jy[k];
	    double z = _vy[k],
		   *_z = _r + _iy[k] * n;
	    for (int i = _nx[j-1]; i < _nx[j]; i++)
		_z[_ix[i]] += _vx[i] * z;
	}
    }
	
#ifdef _TIME_H
    t2 = clock();
    if (R_verbose && *LOGICAL(R_verbose))
	Rprintf("tcrossprod_stm_stm: %.3fs [%.3fs/%.3fs]\n",
		((double) t2 - t0) / CLOCKS_PER_SEC,
		((double) t1 - t0) / CLOCKS_PER_SEC,
		((double) t2 - t1) / CLOCKS_PER_SEC);
#endif
    UNPROTECT(4);
    if (vx != VECTOR_ELT(x, 2))
	UNPROTECT(1);
    if (!isNull(y) &&
	vy != VECTOR_ELT(y, 2))
	UNPROTECT(1);

    return r;
}

// tcrossprod for some triplet matrix and matrix
//
// NOTES 1) tcrossprod does not implement na.rm, so neither do we.
//       2) triplet on triplet does not fit in here.
//       3) if y contains special values we call some bailout
//          function.
//       4) pkgEnv = NULL deactivates the bailout.
//       5) transpose 
//
SEXP tcrossprod_stm_matrix(SEXP x, SEXP R_y, SEXP pkgEnv, SEXP R_verbose,
							  SEXP R_transpose) {
    if (isNull(R_y))
	return tcrossprod_stm_stm(x, R_y, pkgEnv, R_verbose);
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class simple_triplet_matrix");
    if (!isMatrix(R_y))
	error("'y' not of class matrix");

    int n, m;
    SEXP y = R_y;

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
	    r = eval(PROTECT(LCONS(install(".tcrossprod_bailout"),
		     PROTECT( CONS(x,
			      CONS(y, 
			      CONS((R_transpose && *LOGICAL(R_transpose)) ?
				    R_transpose : ScalarLogical(FALSE),
				    R_NilValue)))))), pkgEnv);
	    UNPROTECT(2);
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
	UNPROTECT(2);		/* v, r */
	PROTECT(r);
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




// test validity of list components.
int _valid_ssa(SEXP x) {
    if (LENGTH(x) < 3)
	error("invalid number of components");
    SEXP s = getAttrib(x, R_NamesSymbol);
    int ok = 
	strcmp(CHAR(STRING_ELT(s, 0)), "i") ||
	strcmp(CHAR(STRING_ELT(s, 1)), "v") ||
	strcmp(CHAR(STRING_ELT(s, 2)), "dim") ||
    ((LENGTH(s) > 3) ?
	strcmp(CHAR(STRING_ELT(s, 3)), "dimnames") : 0);
    if (!ok) {
	if (TYPEOF(VECTOR_ELT(x, 0)) != INTSXP ||
	    TYPEOF(VECTOR_ELT(x, 2)) != INTSXP)
	    error("'i, dim' invalid type");
	if (!isVector(VECTOR_ELT(x, 1)))
	    error("'v' not a vector");
	int *xi, *xd, nr, nc;
	s = VECTOR_ELT(x, 0);
	if (!isMatrix(s))
	    error("'i' not a matrix");
	xi = INTEGER(s);
	s  = getAttrib(s, R_DimSymbol);
	nr = INTEGER(s)[0];
	if (nr != LENGTH(VECTOR_ELT(x, 1)))
	    error("'i, v' invalid length");
	nc = INTEGER(s)[1];
	s  = VECTOR_ELT(x, 2);
	if (nc != LENGTH(s))
	    error("'i, dim' invalid length");
	xd = INTEGER(s);
	for (int j = 0; j < nc; j++) {
	    int n = xd[j];
	    if (n > 0) {
		if (n == NA_INTEGER)
		    error("'dim' invalid");
		for (int i = 0; i < nr; i++)
		    if (xi[i] < 1 || xi[i] > n)
			error("i invalid");
	    } else
		if (n < 0)
		    error("'dim' invalid");
		else
		    if (nr > 0)
			error("'dim, i' invalid number of rows");
	    xi += nr;
	}
	if (LENGTH(x) > 3) {
	    s = VECTOR_ELT(x, 3);
	    if (!isNull(s)) {
		if (TYPEOF(s) != VECSXP)
		    error("'dimnames' invalid type");
		if (LENGTH(s) != nc)
		    error("'dimnames' invalid length");
		for (int j = 0; j < nc; j++)
		    if (!isNull(VECTOR_ELT(s, j)) &&
			(LENGTH(VECTOR_ELT(s, j)) != xd[j] ||
		      !isString(VECTOR_ELT(s, j))))
		    error("'dimnames' component invalid length or type");
	    }
	}
    }
    return ok;
}

// wrapper
SEXP __valid_ssa(SEXP x) {
    if (!inherits(x, "simple_sparse_array"))
	return ScalarLogical(FALSE);
    return ScalarLogical(_valid_ssa(x) == FALSE);
}


//
