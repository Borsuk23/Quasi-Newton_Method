#include <math.h>
#include "Equations.h"

integer N = 3;
doublereal x[3] = {1, 1, 1};
doublereal F[3] = { 0, 0, 0 };
doublereal Bk[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
doublereal dk[3] = { 1, 1, 1 };

//doublereal FunctionValue[2];

getFunction(doublereal* FunctionValue, doublereal* x)
{
	FunctionValue[0] = (double)( 4*sin(x[0]) + x[1] + x[2] );
	FunctionValue[1] = (double)( x[1] - 5 );
	FunctionValue[2] = (double)( x[2] + 3 );
}