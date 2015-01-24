#include <math.h>
#include "Equations.h"

integer N = 8;
//doublereal x[3] = {1, 1, 1};
//doublereal F[3] = { 0, 0, 0 };
//doublereal Bk[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//doublereal dk[3] = { 1, 1, 1 };
//doublereal dx[3] = { 0, 0, 0 };

//doublereal FunctionValue[2];


void initData()
{
	int i;
	Bk = (doublereal*)calloc(N*N, sizeof(doublereal));
	//Bk jako macierz jednosktowa
	for (i = 0; i < N*N; i += N + 1)
	{
		Bk[i] = 1;
	}
	F = (doublereal*)calloc(N, sizeof(doublereal));
	dk = (doublereal*)calloc(N, sizeof(doublereal));
	//dk={1,1,1,..}
	for (i = 0; i < N; i++)
	{
		dk[i] = 1;
	}
	x = (doublereal*)calloc(N, sizeof(doublereal));
	dx = (doublereal*)calloc(N, sizeof(doublereal));
	dF = (doublereal*)calloc(N, sizeof(doublereal));
	ipiv = (integer*)calloc(N, sizeof(integer));
	licznik = (doublereal*)calloc(N, sizeof(doublereal));
	mianownik[0] = 1;
}

void clearData()
{
	free(Bk);
	free(F);
	free(dk);
	free(x);
	free(dx);
	free(dF);
	free(ipiv);
	free(licznik);
}


getFunction(doublereal* FunctionValue, doublereal* x)
{
	//FunctionValue[0] = (double)( 4*sin(x[0]) + x[1] + x[2] );
	//FunctionValue[1] = (double)( x[1] - 5 );
	//FunctionValue[2] = (double)( x[2] + 3 );

	//FunctionValue[0] = (double)(pow(x[0], 2) - pow(x[1], 2) + x[2] );
	//FunctionValue[1] = (double)(pow(x[0], 2) - pow(x[1], 2) + pow(x[2], 2) - 6 );
	//FunctionValue[2] = (double)(x[0] + x[1] - 2 * pow(x[2], 2) + 15);

	//FunctionValue[0] = (double)(x[0] - x[1] + x[2] - 4);
	//FunctionValue[1] = (double)(x[0] - x[2] + 3);
	//FunctionValue[2] = (double)(x[2] - x[0] - x[1] );

	//Mikolaj
	FunctionValue[0] = (double)(pow(x[0], 5) + pow(x[1], 3) + pow(x[2], 2) + pow(x[3], 2) + pow(x[4], 2) + pow(x[5], 2) - pow(x[6], 2) - pow(x[7], 2) - 1177);
	FunctionValue[1] = (double)(pow(x[0], 1) + pow(x[1], 1) + pow(x[2], 1) + pow(x[3], 1) + pow(x[4], 1) + pow(x[5], 1) - pow(x[6], 1) - pow(x[7], 1) - 45);
	FunctionValue[2] = (double)(pow(x[1], 1) - pow(x[2], 1) + pow(x[3], 1) + pow(x[4], 1) + 1);
	FunctionValue[3] = (double)(pow(x[2], 1) + pow(x[3], 1) - 22);
	FunctionValue[4] = (double)(pow(x[3], 1) + pow(x[4], 1) + pow(x[5], 1) - 43);
	FunctionValue[5] = (double)(pow(x[4], 1) + pow(x[5], 1) - 35);
	FunctionValue[6] = (double)(pow(x[4], 1) + pow(x[5], 1) + pow(x[6], 1) - 45);
	FunctionValue[7] = (double)(pow(x[6], 1) + pow(x[7], 1) - 15);
	
	//FunctionValue[0] = (double)(cos(x[0]) + x[1]);
	//FunctionValue[1] = (double)(x[0] * x[1]);

	//FunctionValue[0] = (double)(x[0] - 8);
	//FunctionValue[2] = (double)(x[2] + 12);
	//FunctionValue[1] = (double)(x[1] - 12 * x[2] - 4);

	//FunctionValue[0] = (double)(2 * x[0] - x[1] + x[2] + 4);
	//FunctionValue[1] = (double)(8 * x[0] + 2 * x[1] - 5 * x[2] + 10);
	//FunctionValue[2] = (double)(4 * x[0] + x[1] + x[2] - 2);

	/*FunctionValue[0] = (double)( x[0] + x[1] + x[2] );
	FunctionValue[1] = (double)( x[4] - x[3] + 10);
	FunctionValue[2] = (double)( 6 - x[3] - x[2] );
	FunctionValue[3] = (double)( 24 + 6*x[7] );
	FunctionValue[4] = (double)( x[5] + x[0] );
	FunctionValue[5] = (double)( 2 - x[0]);
	FunctionValue[6] = (double)( x[2] - 1 );
	FunctionValue[7] = (double)( 4*x[6]);*/
}