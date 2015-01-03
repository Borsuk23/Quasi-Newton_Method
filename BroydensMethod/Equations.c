#include <math.h>
#include "Equations.h"

integer N = 3;
//doublereal x[3] = {1, 1, 1};
doublereal F[3] = { 0, 0, 0 };
doublereal Bk[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
doublereal dk[3] = { 1, 1, 1 };

//doublereal FunctionValue[2];

getFunction(doublereal* FunctionValue, doublereal* x)
{
	//FunctionValue[0] = (double)( 4*sin(x[0]) + x[1] + x[2] );
	//FunctionValue[1] = (double)( x[1] - 5 );
	//FunctionValue[2] = (double)( x[2] + 3 );

	//FunctionValue[0] = (double)(pow(x[0], 1) - pow(x[1], 1) + x[2] - 8);
	//FunctionValue[1] = (double)(pow(x[0], 1) - pow(x[1], 1) + pow(x[2], 1) + 12);
	//FunctionValue[2] = (double)(x[0] + pow(x[1], 1) - 12 * pow(x[2], 1) - 4);

	FunctionValue[0] = (double)(x[0] - x[1] + x[2] - 4);
	FunctionValue[1] = (double)(x[0] - x[2] + 3);
	FunctionValue[2] = (double)(x[2] - x[0] - x[1] );

	//FunctionValue[0] = (double)(x[0] - 8);
	//FunctionValue[2] = (double)(x[2] + 12);
	//FunctionValue[1] = (double)(x[1] - 12 * x[2] - 4);

	//FunctionValue[0] = (double)(2 * x[0] - x[1] + x[2] + 4);
	//FunctionValue[1] = (double)(8 * x[0] + 2 * x[1] - 5 * x[2] + 10);
	//FunctionValue[2] = (double)(4 * x[0] + x[1] + x[2] - 2);


}