#include <math.h>
#include "Equations.h"

integer N = 2;
doublereal x[2] = {1, 1};
doublereal F[2] = { 0, 0 };
doublereal Bk[4] = { 1, 0, 0, 1 };
doublereal dk[2] = { 1, 1 };

//doublereal FunctionValue[2];

getFunction(doublereal* FunctionValue, doublereal* x)
{
	FunctionValue[0] = (double)( cos(x[0]) + x[1] );
	FunctionValue[1] = (double)( x[0] * x[1] );
}