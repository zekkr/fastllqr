#include <R.h>
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern void F77_NAME(llqr_ppro_fortran)(void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,
                                        void *, void *, void *);

extern void F77_NAME(tvcqr_seq_ppro_fortran)(void *, void *, void *, void *, void *,
                                             void *, void *, void *, void *, void *,
                                             void *, void *, void *, void *, void *,
                                             void *, void *, void *, void *, void *,
                                             void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"llqr_ppro_fortran", (DL_FUNC) &F77_NAME(llqr_ppro_fortran), 23},
    {"tvcqr_seq_ppro_fortran", (DL_FUNC) &F77_NAME(tvcqr_seq_ppro_fortran), 21},
    {NULL, NULL, 0}
};

void attribute_visible R_init_fastllqr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
