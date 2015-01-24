#include <math.h>
#include "Equations.h"

integer N = 8;



void initData()
{
	int i;
	Bk = (doublereal*)calloc(N*N, sizeof(doublereal));
	Bkpom = (doublereal*)calloc(N*N, sizeof(doublereal));
	dxpom = (doublereal*)calloc(N*N, sizeof(doublereal));
	eye = (doublereal*)calloc(N*N, sizeof(doublereal));
	//Bk jako macierz jednosktowa
	for (i = 0; i < N*N; i += N + 1)
	{
		Bk[i] = 1;
		Bkpom[i] = 1;
		dxpom[i] = 1;
		eye[i] = 1;
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
	work = (doublereal*)calloc(2 * N, sizeof(doublereal));
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
	//FunctionValue[0] = (double)(x[0] + pow(x[1], 2) - 5.5*x[1] + pow(x[2], 2) - 11 * x[2] / 3 + 35 / 6);
	//FunctionValue[1] = (double)(-pow(x[0], 2) + 10.5*x[0] - pow(x[1], 2) + 19 * x[1] / 3 + pow(x[2], 2) - 3.75*x[2] - 34);
	//FunctionValue[2] = (double)(pow(x[0], 3) - 6 * pow(x[0], 2) - 44 * x[0] / 3 + 0.25*x[1] + 0.2*x[2] + 5831 / 60);
	
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

void initBroyden(doublereal* function, doublereal* x, float dokladnosc)
{
	int i = 0;
	int j = 0;
	int k = 0;
	getFunction(function, x);
	doublereal* dx = (doublereal*)calloc(N, sizeof(doublereal));
	doublereal* deltaFunction = (doublereal*)calloc(N, sizeof(doublereal));
	for (i = 0; i < N; i++)
	{
		for (k = 0; k < N; k++)
		{
			if (k == i)
			{
				dx[k] = x[k] + dokladnosc;
			}
			else
			{
				dx[k] = x[k];
			}
		}
		getFunction(deltaFunction, dx);
		for (j = i*N; j < N + i*(N); j++)
		{
			Bk[j] = (deltaFunction[j%N] - function[j%N]) / dokladnosc;
		}
	}
	printf("Macierz Bk:\n");
	for (i = 0; i < N; i++)
	{
		printf("      [");
		for (integer j = i; j < (N*N); j += N)
		{
			printf(" %lf ", Bk[j]);
		}
		printf("]\n");
	}
	printf("\n");
	getFunction(F, x);
	free(dx);
	free(deltaFunction);
}
