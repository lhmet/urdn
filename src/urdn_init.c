#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(fluxdir)(void *, void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *, void *,
                     void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"fluxdir", (DL_FUNC) &F77_NAME(fluxdir), 15},
  {NULL, NULL, 0}
};

void R_init_urdn(DllInfo *dll)
{
  /*                 dll,          .C,            .Call,           .Fortran,           .External */
  /*                 dll, R_CMethodDef, R_CallMethodDef, R_FortranMethodDef, R_ExternalMethodDef */
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
