#include "f2c.h"
#include "clapack.h"

integer N;
doublereal x[2];
doublereal F[2];
doublereal Bk[4];
doublereal dk[2];
doublereal dx[2];
doublereal dF[2];
integer ipiv[2];

getFunction(doublereal*, doublereal*);