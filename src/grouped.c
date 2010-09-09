#include <R.h>
#include <Rdefines.h>
#include <time.h>

extern int _valid_stm(SEXP x);

// ceeboo 2010/8
//
// sum (collapse) the rows of x into the column groups 
// defined in index.
//
SEXP _row_tsums(SEXP x, SEXP R_index, SEXP R_na_rm, SEXP R_verbose) {
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class 'simple_triplet_matrix'");
    if (!inherits(R_index, "factor"))
	error("'index' not of class 'factor'");

    int *p, *q, k, n, m, l;
    SEXP _v, _i, _j, __v, r, s;

    if (LENGTH(R_index) != INTEGER(VECTOR_ELT(x, 4))[0])
	error("'index' invalid length");
    //m = ;
    m = INTEGER(VECTOR_ELT(x, 3))[0];

#ifdef _TIME_H
    // code section times
    clock_t t2, t1, t0 = clock();
#endif
    _i = VECTOR_ELT(x, 0);

    p = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(_i))));
    q = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(_i))));

    _j = VECTOR_ELT(x, 1);

    // NOTE tapply ignores NA(s) in argument INDEX.
    //      So do we with a warning.
    //      
    // FIXME bailout if k * m > 2147483647.
    l = 0;
    for (int i = 0; i < LENGTH(_i); i++) {
	k = INTEGER(R_index)[INTEGER(_j)[i] - 1];
	if (k == NA_INTEGER)
	    continue;
	k = INTEGER(_i)[i] - 1 + (k - 1) * m;
	p[l] = k;
	q[l] = i;
	l++;
    }
    if (l < LENGTH(_i))
	warning("NA(s) in 'index'");
    else
	for (int i = 0; i < LENGTH(R_index); i++)
	    if (INTEGER(R_index)[i] == NA_INTEGER) {
		warning("NA(s) in 'index'");
		break;
	    }

    if (l)
	R_qsort_int_I(p, q, 1, l);

    n =  0;
    k = -1;
    for (int i = 0; i < l; i++)
	if (p[i] > k) {
	    k = p[i];
	    n++;
	}

    // NOTE use REALSXP to avoid overflows.
    r = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(r, 0, ( _i = allocVector(INTSXP,  n)));
    SET_VECTOR_ELT(r, 1, ( _j = allocVector(INTSXP,  n)));
    SET_VECTOR_ELT(r, 2, (__v = allocVector(REALSXP, n)));
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

    n = -1;
    k = -1;
    switch (TYPEOF(_v)) {
	case LGLSXP:
	case INTSXP:
	    for (int i = 0; i < l; i++) {
		if (p[i] > k) {
		    k = p[i];
		    n++;
		    INTEGER(_i)[n] =  k % m + 1;
		    INTEGER(_j)[n] =  k / m + 1;
		    REAL(__v)[n]   =  0;
		}
		int z = INTEGER(_v)[q[i]];
		if (z != NA_INTEGER)
		    REAL(__v)[n] += (double) z;
		else
		    if (!LOGICAL(R_na_rm)[0])
			REAL(__v)[n] += NA_REAL;
	    }
	    break;
	case REALSXP:
	    for (int i = 0; i < l; i++) {
		if (p[i] > k) {
		    k = p[i];
		    n++;
		    // NOTE we expect that retrieving the i
		    //      index value is less efficient
		    INTEGER(_i)[n] =  k % m + 1;
		    INTEGER(_j)[n] =  k / m + 1;
		    REAL(__v)[n]   =  0;
		}
		double z = REAL(_v)[q[i]];
		// NOTE with na.rm == TRUE the sum may be zero
		if (!ISNAN(z) || !LOGICAL(R_na_rm)[0])
		    REAL(__v)[n] += z;
	    }
	    break;
	default:
	    error("type of 'v' not supported");
    }
#ifdef _TIME_H
    t2 = clock();
    if (R_verbose && *LOGICAL(R_verbose))
        Rprintf("_row_tsums: %.3fs [%.3fs/%.3fs]\n", 
                ((double) t2 - t0) / CLOCKS_PER_SEC,
                ((double) t1 - t0) / CLOCKS_PER_SEC,
                ((double) t2 - t1) / CLOCKS_PER_SEC);
#endif
    UNPROTECT(3);

    return r;
}
