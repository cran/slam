
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP __valid_stm(SEXP x);
extern SEXP __valid_ssa(SEXP x);
extern SEXP __valid_v(SEXP x);
extern SEXP _split_col(SEXP x);
extern SEXP _all_row(SEXP x, SEXP _na_rm);
extern SEXP _part_index(SEXP x);
extern SEXP _vector_index(SEXP d, SEXP x);
extern SEXP _ini_array(SEXP d, SEXP p, SEXP v, SEXP s);
extern SEXP _match_matrix(SEXP x, SEXP y, SEXP _nm);
extern SEXP _unattr(SEXP x);
extern SEXP _sums_stm(SEXP x, SEXP R_dim, SEXP R_na_rm);
extern SEXP _row_tsums(SEXP x, SEXP R_index, SEXP R_na_rm, SEXP R_reduce, 
		       SEXP R_verbose);
extern SEXP tcrossprod_stm_stm(SEXP x, SEXP y, SEXP pkgEnv, SEXP R_verbose);
extern SEXP tcrossprod_stm_matrix(SEXP x, SEXP R_y, SEXP pkgEnv, 
				  SEXP R_verbose, SEXP R_transpose);
extern SEXP _col_apply_stm(SEXP a);


static const R_CallMethodDef CallEntries[] = {
    {"R__valid_stm",		(DL_FUNC) __valid_stm,		 1},
    {"R__valid_ssa",		(DL_FUNC) __valid_ssa,		 1},
    {"R__valid_v",		(DL_FUNC) __valid_v,		 1},
    {"R_split_col",		(DL_FUNC) _split_col,		 1},
    {"R_all_row",		(DL_FUNC) _all_row,		 2},
    {"R_part_index",		(DL_FUNC) _part_index,		 1},
    {"R_vector_index",		(DL_FUNC) _vector_index,	 2},
    {"R_ini_array",		(DL_FUNC) _ini_array,		 4},
    {"R_match_matrix",		(DL_FUNC) _match_matrix,	 3},
    {"R_unattr",		(DL_FUNC) _unattr,		 1},
    {"R_sums_stm",		(DL_FUNC) _sums_stm,		 3},
    {"R_row_tsums",		(DL_FUNC) _row_tsums,		 5},
    {"R_tcrossprod_stm_matrix", (DL_FUNC) tcrossprod_stm_matrix, 5},
    {"R_tcrossprod_stm_stm",	(DL_FUNC) tcrossprod_stm_stm,	 4},
    {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExternalEntries[] = {
    {"R_col_apply_stm",		(DL_FUNC) _col_apply_stm,	-1},
    {NULL, NULL, 0}
};

void R_init_slam(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, ExternalEntries);
    R_useDynamicSymbols(dll, FALSE);
}

