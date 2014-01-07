
#include <R.h>
#include <Rdefines.h>
#include <time.h>

extern int _valid_stm(SEXP x);

// (C) ceeboo 2013/12
//
// Wrapper for simple triplet matrix which runs in constant 
// memory.
SEXP _col_apply_stm(SEXP a) {
    a = CDR(a);
    if (length(a) < 2)
	error("invalid number of arguments");

    SEXP x = CAR(a);
    if (!inherits(x, "simple_triplet_matrix") || _valid_stm(x))
	error("'x' not of class 'simple_triplet_matrix'");
    
    if (!isFunction(CADR(a)))
	error("invalid function parameter");

    int n, m, *_ix, *_jx, *_nx, *_px;
    SEXP vx, z, r;

    vx = VECTOR_ELT(x, 2);

    n = INTEGER(VECTOR_ELT(x, 3))[0];
    m = INTEGER(VECTOR_ELT(x, 4))[0];

    z = PROTECT(allocVector(TYPEOF(vx), n));
    a =	PROTECT(LCONS(CADR(a), LCONS(z, CDDR(a))));

    switch(TYPEOF(vx)) {
        case LGLSXP:
        case INTSXP:
            memset(INTEGER(z), 0, sizeof(int) * n);
            break;
        case REALSXP:
            memset(REAL(z), 0, sizeof(double) * n);
            break;
        case RAWSXP:
            memset(RAW(z), 0, sizeof(char) * n);
            break;
        case CPLXSXP:
            memset(COMPLEX(z),	0, sizeof(Rcomplex) * n);
            break;
        case EXPRSXP:
        case VECSXP:
            for (int i = 0; i < n; i++)
                SET_VECTOR_ELT(z, i, R_NilValue);
            break;
        case STRSXP:
            for (int i = 0; i < n; i++)
                SET_STRING_ELT(z, i, R_BlankString);
            break;
        default:
            error("type of 'v' not supported");
    }

    // Map blocks of equal column indexes

    _jx = INTEGER(VECTOR_ELT(x, 1));	// column indexes

    _nx = INTEGER(PROTECT(allocVector(INTSXP, m + 1)));
    memset(_nx, 0, sizeof(int) * (m + 1));
    for (int k = 0; k < LENGTH(vx); k++)
	_nx[_jx[k]]++;
    for (int k = 1; k < m + 1; k++)
	_nx[k] += _nx[k-1];

    _px = INTEGER(PROTECT(allocVector(INTSXP, LENGTH(vx))));
    _nx -= 1;				// one-based R indexing
    for (int k = 0; k < LENGTH(vx); k++) {
	_px[_nx[_jx[k]]] = k;
	_nx[_jx[k]]++;	    
    }
    // Reset
    _nx += 1;
    for (int k = m; k > 0; k--)
	_nx[k] = _nx[k-1];
    _nx[0] = 0;

    _ix = INTEGER(VECTOR_ELT(x, 0));    // row indexes

    r = PROTECT(allocVector(VECSXP, m));

    int f, fl;
    f = fl = _nx[0];
    for (int i = 1; i < m + 1; i++) {
	int l = _nx[i];
	// (Re)set values
	switch(TYPEOF(vx)) {
	    case LGLSXP:
	    case INTSXP:
		for (int k = fl; k < f; k++)
		    INTEGER(z)[_px[k]] = 0;
		for (int k = f; k < l; k++) {
		    int p = _px[k],
			i = _ix[p] - 1;
		    INTEGER(z)[i] = INTEGER(vx)[p];
		    _px[k] = i;
		}
		break;
	    case REALSXP:
		for (int k = fl; k < f; k++)
		    REAL(z)[_px[k]] = 0.0;
		for (int k = f; k < l; k++) {
		    int p = _px[k],
			i = _ix[p] - 1;
		    REAL(z)[i] = REAL(vx)[p];
		    _px[k] = i;
		}
		break;
	    case RAWSXP:
		for (int k = fl; k < f; k++)
		    RAW(z)[_px[k]] = (char) 0;
		for (int k = f; k < l; k++) {
		    int p = _px[k],
			i = _ix[p] - 1;
		    RAW(z)[i] = RAW(vx)[p];
		    _px[k] = i;
		}
		break;
	    case CPLXSXP:
		for (int k = fl; k < f; k++) {
		    static Rcomplex c;
		    COMPLEX(z)[_px[k]] = c;
		}
		for (int k = f; k < l; k++) {
		    int p = _px[k],
			i = _ix[p] - 1;
		    COMPLEX(z)[i] = COMPLEX(vx)[p];
		    _px[k] = i;
		}
		break;
	    case EXPRSXP:
	    case VECSXP:
		for (int k = fl; k < f; k++)
		    SET_VECTOR_ELT(z, _px[k], R_NilValue);
		for (int k = f; k < l; k++) {
		    int p = _px[k],
			i = _ix[p] - 1;
		    SET_VECTOR_ELT(z, i, VECTOR_ELT(vx, p));
		    _px[k] = i;
		}
		break;
	    case STRSXP:
		for (int k = fl; k < f; k++)
		    SET_STRING_ELT(z, _px[k], R_BlankString);
		for (int k = f; k < l; k++) {
		    int p = _px[k],
			i = _ix[p] - 1;
		    SET_STRING_ELT(z, i, STRING_ELT(vx, p));
		    _px[k] = i;
		}
		break;
	    default:
		error("type of 'v' not supported");
	}
	SEXP s = eval(a, R_GlobalEnv);
	if (s == z)			// identity, print, ...
	    SET_VECTOR_ELT(r, i - 1, duplicate(s));
	else
	    SET_VECTOR_ELT(r, i - 1, s);
	fl = f;
	f  = l;
    }

    UNPROTECT(5);
    return r;
}

