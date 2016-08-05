#include <R.h>
#include <Rdefines.h>
#include <time.h>

extern int _valid_stm(SEXP x);

// ceeboo 2010/8+10, 2016/6
//
// sum (collapse) the rows of x into the column groups 
// defined in index.
//
SEXP _row_tsums(SEXP x, SEXP R_index, SEXP R_na_rm, SEXP R_reduce, 
		SEXP R_verbose) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class 'simple_triplet_matrix'");
    if (!inherits(R_index, "factor"))
	error("'index' not of class 'factor'");

    int *p, *q, k, n, m, f, l;
    SEXP _v, _i, _j, __i, __v, r, s;

    if (LENGTH(R_index) != INTEGER(VECTOR_ELT(x, 4))[0])
	error("'index' invalid length");
    if (TYPEOF(R_na_rm) != LGLSXP)
	error("'na_rm' not logical");
    if (!LENGTH(R_na_rm))
	error("'na_rm' invalid length");
    int na_rm = LOGICAL(R_na_rm)[0] == TRUE;

    if (TYPEOF(R_reduce) != LGLSXP)
	error("'reduce' not logical");
    if (!LENGTH(R_reduce))
	error("'reduce' invalid length");

#ifdef _TIME_H
    // code section times
    clock_t t2, t1, t0 = clock();
#endif
    _i = VECTOR_ELT(x, 0);

    p = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(_i))));
    q = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(_i))));

    // sort by row indexes
    for (int i = 0; i < LENGTH(_i); i++) {
	p[i] = INTEGER(_i)[i];
	q[i] = i;
    }
    if (LENGTH(_i))
	R_qsort_int_I(p, q, 1, LENGTH(_i));

    // sort row blocks by column indexes
    //
    // NOTE we change the sign with each block
    //	    to ensure a change in key.
    //
    _j = VECTOR_ELT(x, 1);

    f = 0;
    l = 0;
    n = 0;
    m = 0;
    for (int i = 0; i < LENGTH(_i); i++) {
	k = INTEGER(R_index)[INTEGER(_j)[q[i]] - 1];
	if (k == NA_INTEGER)
	    continue;
	if (n != p[i]) {
	    n  = p[i];
	    if (f < l)
		R_qsort_int_I(p, q, f, l);
	    f = l + 1;
	    m = (m) ? 0 : 1;
	}
	p[l] = (m) ? k : -k;
	q[l] = q[i];
	l++;
    }
    if (l) {
	R_qsort_int_I(p, q, f, l);
	// FIXME this may be time-consuming.
	if (l < LENGTH(_i))
	    warning("NA(s) in 'index'");
	else
	    for (int i = 0; i < LENGTH(R_index); i++)
		if (INTEGER(R_index)[i] == NA_INTEGER) {
		    warning("NA(s) in 'index'");
	    	    break;
		}
    }

    // count
    n = 0;
    k = 0;
    for (int i = 0; i < l; i++)
	if (k != p[i]) {
	    k  = p[i];
	    n++;
	}

    r = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(r, 0, (__i = allocVector(INTSXP,  n)));
    SET_VECTOR_ELT(r, 1, ( _j = allocVector(INTSXP,  n)));
    SET_VECTOR_ELT(r, 3, VECTOR_ELT(x, 3));
    SET_VECTOR_ELT(r, 4, 
	ScalarInteger(LENGTH(getAttrib(R_index, R_LevelsSymbol))));

    SET_VECTOR_ELT(r, 5, (s = allocVector(VECSXP, 2)));
    SET_VECTOR_ELT(s, 0, R_NilValue);
    SET_VECTOR_ELT(s, 1, getAttrib(R_index, R_LevelsSymbol));
    if (LENGTH(x) > 5) {
	SEXP t = VECTOR_ELT(x, 5);
	if (!isNull(t)) {
	    SET_VECTOR_ELT(s, 0, VECTOR_ELT(t, 0));
	    if (!isNull((t = getAttrib(t, R_NamesSymbol)))) 
		setAttrib(s, R_NamesSymbol, t);
	}
	setAttrib(r, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    }
    else {
	setAttrib(r, R_NamesSymbol, (s = allocVector(STRSXP, 6)));
	SEXP t = getAttrib(x, R_NamesSymbol);
	for (int i = 0; i < 5; i++)
	    SET_STRING_ELT(s, i, STRING_ELT(t, i));
	SET_STRING_ELT(s, 5, mkString("dimnames"));
    }
    setAttrib(r, R_ClassSymbol, getAttrib(x, R_ClassSymbol));
#ifdef _TIME_H
    t1 = clock();
#endif
    _v = VECTOR_ELT(x, 2);

    switch (TYPEOF(_v)) {
	case LGLSXP:
	case INTSXP:
	    {
		// NOTE use REALSXP to avoid overflows.
		SET_VECTOR_ELT(r, 2, (__v = allocVector(REALSXP, n)));
		double *_z = NULL;

		n = 0;
		k = 0;
		for (int i = 0; i < l; i++) {
		    if (k != p[i]) {
			k  = p[i];
			INTEGER(__i)[n] = INTEGER(_i)[q[i]];
			INTEGER( _j)[n] = (k > 0) ? k : -k;
			 _z = REAL(__v) + n;    
			*_z = 0;
			n++;
		    }
		    int z = INTEGER(_v)[q[i]];
		    if (z != NA_INTEGER)
			*_z += (double) z;
		    else
			if (!na_rm)
			    *_z = NA_REAL;
		}
	    }
	    break;
	case REALSXP:
	    {
		SET_VECTOR_ELT(r, 2, (__v = allocVector(REALSXP, n)));
		double *_z = NULL;

		n = 0;
		k = 0;
		for (int i = 0; i < l; i++) {
		    if (k != p[i]) {
			k  = p[i];
			INTEGER(__i)[n] = INTEGER(_i)[q[i]];
			INTEGER( _j)[n] = (k > 0) ? k : -k;
			 _z = REAL(__v) + n;    
			*_z = 0;
			n++;
		    }
		    double z = REAL(_v)[q[i]];
		    if (!na_rm || !ISNAN(z))
			*_z += z;
		}
	    }
	    break;
	case CPLXSXP:
	    {
		SET_VECTOR_ELT(r, 2, (__v = allocVector(CPLXSXP, n)));
		Rcomplex *_z = NULL;

		n = 0;
		k = 0;
		for (int i = 0; i < l; i++) {
		    if (k != p[i]) {
			k  = p[i];
			INTEGER(__i)[n] = INTEGER(_i)[q[i]];
			INTEGER( _j)[n] = (k > 0) ? k : -k;
			_z = COMPLEX(__v) + n;    
			_z->r = 0;
			_z->i = 0;
			n++;
		    }
		    Rcomplex *z = COMPLEX(_v) + q[i];
		    if (!na_rm || (!ISNAN(z->r) && !ISNAN(z->i))) {
			_z->r += z->r;
			_z->i += z->i;
		    }
		}
	    }
	    break;
	default:
	    error("type of 'v' invalid");
    }

    // remove zeros
    if (*LOGICAL(R_reduce)) {
	k = n;
	n = 0;
	if (TYPEOF(__v) == CPLXSXP)
	    for (int i = 0; i < k; i++) {
		if (COMPLEX(__v)[i].r == 0.0 && 
		    COMPLEX(__v)[i].i == 0.0)
		    continue;
		if (i > n) {
		    INTEGER(__i)[n] = INTEGER(__i)[i];
		    INTEGER( _j)[n] = INTEGER( _j)[i];
		    COMPLEX(__v)[n] = COMPLEX(__v)[i];
		}
		n++;
	    }
	else
	    for (int i = 0; i < k; i++) {
		if (REAL(__v)[i] == 0.0)
		    continue;
		if (i > n) {
		    INTEGER(__i)[n] = INTEGER(__i)[i];
		    INTEGER( _j)[n] = INTEGER( _j)[i];
		       REAL(__v)[n] =    REAL(__v)[i];
		}
		n++;
	    }
	if (n < k) {
	    SETLENGTH(__i, n);
	    SETLENGTH( _j, n);
	    SETLENGTH(__v, n);
	}
    }
#ifdef _TIME_H
    t2 = clock();
    if (R_verbose && *LOGICAL(R_verbose)) {
	if (*LOGICAL(R_reduce))
	    Rprintf("_row_tsums: reduced %i (%i) zeros\n", k - n, n);
        Rprintf("_row_tsums: %.3fs [%.3fs/%.3fs]\n", 
                ((double) t2 - t0) / CLOCKS_PER_SEC,
                ((double) t1 - t0) / CLOCKS_PER_SEC,
                ((double) t2 - t1) / CLOCKS_PER_SEC);
    }
#endif
    UNPROTECT(3);

    return r;
}

