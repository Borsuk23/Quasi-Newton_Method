#include "f2c.h"
#include "clapack.h"

integer N;
doublereal x[3];
doublereal F[3];
doublereal Bk[9];
doublereal dk[3];
doublereal dx[3];
doublereal dF[3];
integer ipiv[3];

getFunction(doublereal*, doublereal*);