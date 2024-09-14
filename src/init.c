#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
extern SEXP C_ComputeLD(SEXP, SEXP);
extern SEXP C_Diff(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_IsAscending(SEXP, SEXP, SEXP);
extern SEXP C_MakeGrid(SEXP, SEXP, SEXP);
extern SEXP C_NullD(SEXP, SEXP);
extern SEXP C_SbarBlocks(SEXP, SEXP, SEXP);
extern SEXP C_SbarLTB(SEXP, SEXP);
static const R_CallMethodDef CallEntries[] = {
    {"C_ComputeLD", (DL_FUNC) &C_ComputeLD, 2},
    {"C_Diff", (DL_FUNC) &C_Diff, 4},
    {"C_IsAscending", (DL_FUNC) &C_IsAscending, 3},
    {"C_MakeGrid", (DL_FUNC) &C_MakeGrid, 3},
    {"C_NullD", (DL_FUNC) &C_NullD, 2},
    {"C_SbarBlocks", (DL_FUNC) &C_SbarBlocks, 3},
    {"C_SbarLTB", (DL_FUNC) &C_SbarLTB, 2},
    {NULL, NULL, 0}
};
void R_init_gps_mgcv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
