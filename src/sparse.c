#include <R.h>
#include <Rdefines.h>

// ceeboo 2009/5
//


// row or column sums of some triplet matrix
//
SEXP _sums_stm(SEXP x, SEXP R_dim, SEXP R_na_rm) {
    if (!inherits(x, "stm"))
	error("'x' not of class 'stm'");
    if (TYPEOF(R_dim) != INTSXP)
	error("'dim' not of type integer");
    if (TYPEOF(R_na_rm) != LGLSXP)
	error("'na.rm' not of type logical");

    int n, *i = NULL;
    
    switch ((n = *INTEGER(R_dim))) {
	case 1:
	    i = INTEGER(getAttrib(x, install("i")));
	    break;
	case 2:
	    i = INTEGER(getAttrib(x, install("j")));
	    break;
	default:
	    error("'dim' invalid");
    }
    n = INTEGER(getAttrib(x, install("Dim")))[n-1];

    // for the type of the return argument see the behavior
    // of rowSums and colSums for matrix.
    SEXP r = PROTECT(allocVector(REALSXP, n));

    memset(REAL(r), 0, sizeof(double) * n);

    switch (TYPEOF(x)) {
	case LGLSXP:
	case INTSXP:
	    if (*LOGICAL(R_na_rm)) {
		int v;
		for (int k = 0; k < LENGTH(x); k++)
		    if ((v = INTEGER(x)[k]) == NA_INTEGER)
			continue;
		    else 
			REAL(r)[i[k]-1] += (double) v;
	    } else {
		int v;
		for (int k = 0; k < LENGTH(x); k++)
		    REAL(r)[i[k]-1] +=
		    // map NA
			 ((v = INTEGER(x)[k]) == NA_INTEGER) ? NA_REAL : v;
	    }
	    break;
	case REALSXP:
	    if (*LOGICAL(R_na_rm)) {
		double v;
		for (int k = 0; k < LENGTH(x); k++)
		    if (ISNAN((v = REAL(x)[k])))
			continue;
		    else
			REAL(r)[i[k]-1] += v;
	    } else
		for (int k = 0; k < LENGTH(x); k++)
		    REAL(r)[i[k]-1] += REAL(x)[k];
	    break;
	default:
	    error("type of 'x' not supported");
    }

    SEXP d = getAttrib(x, install("Dimnames"));
    if (!isNull(d))
	setAttrib(r, install("names"), VECTOR_ELT(d, *INTEGER(R_dim)-1));

    UNPROTECT(1);

    return r;
}

// tcrossprod for some triplet matrix and matrix
//
// NOTES 1) tcrossprod does not implement na.rm, so neither do we.
//       2) triplet on triplet does not fit in here.
//       3) if y = NULL or contains special values we call some
//          bailout function.
//
SEXP tcrossprod_stm_matrix(SEXP x, SEXP y, SEXP pkgEnv) {
    if (!inherits(x, "stm"))
	error("'x' not of class stm");
    if (isNull(y))
	goto bailout;
    if (!isMatrix(y))
	error("'y' not of class matrix");

    int i, j, k, l, n, m;

    k = INTEGER(getAttrib(x, install("Dim")))[1];

    if (k != INTEGER(getAttrib(y, R_DimSymbol))[1])
	error("the number of columns of 'x' and 'y' do not conform");

    n = INTEGER(getAttrib(x, install("Dim")))[0];
    m = INTEGER(getAttrib(y, R_DimSymbol))[0];

    SEXP r = PROTECT(allocMatrix(REALSXP, n, m));

    memset(REAL(r), 0, sizeof(double) * n * m);

    int *xi, *xj;

    xi = INTEGER(getAttrib(x, install("i")));
    xj = INTEGER(getAttrib(x, install("j")));

    // coercing is in general not storage efficient, and therefore
    // assumes that y is not too large. on the other hand, as the
    // entries of y could be accessed multiple times, casting would
    // not be runtime efficient. if memory footprint is of concern
    // then the program flow should be further switch(ed).
    if (TYPEOF(y) != REALSXP)
	y = coerceVector(y, REALSXP);

    // check for special values
    for (k = 0; k < LENGTH(y); k++)
	if (!R_FINITE(REAL(y)[k])) {
	    UNPROTECT(1);
bailout:
	    return eval(LCONS(install(".tcrossprod.bailout"),
			LCONS(x,
			LCONS(y, R_NilValue))), pkgEnv);
	}

    switch (TYPEOF(x)) {
	case LGLSXP:
	case INTSXP:
	    for (k = 0; k < LENGTH(x); k++) {
		i = xi[k] - 1;
		j = xj[k] - 1;

		double z = (double) INTEGER(x)[k];

		for (l = 0; l < m; l++)
		    REAL(r)[i + l * n] += z * REAL(y)[l + j * m];
	    }
	    break;
	case REALSXP:
	    for (k = 0; k < LENGTH(x); k++) {
		i = xi[k] - 1;
		j = xj[k] - 1;

		double z = REAL(x)[k];

		for (l = 0; l < m; l++)
		    REAL(r)[i + l * n] += z * REAL(y)[l + j * m];
	    }
	    break;
	default:
	    error("type of 'x' not supported");
    }

    // set dimnames and names of dimnames.
    SEXP dn = getAttrib(x, install("Dimnames"));

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

