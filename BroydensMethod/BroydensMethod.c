#include < stdio.h>
#include "f2c.h"
#include "clapack.h"
#include "Equations.h"
#include "SaveToFile.h"

int main(void)
{
	FILE * fp;
    integer info;
	float dokladnosc;
	integer i;
	doublereal alpha = -1;
	doublereal zero = 0;
	integer nrhs = 1;
	doublereal xAlpha = 1;
	integer lda = N;
	integer ldb = N;
	doublereal AlphaK = 1;
	integer counter = 0;
	integer typ = 1;
	char typeN[1] = { 'N' }; //macierz normalna
	char typeT[1] = { 'T' }; //macierz transponowana
	

	initData();

	//pobor dokladnosci i punktu startowego z pliku
	printf("Wczytywanie parametrow z pliku..\n");
	fp = fopen("dane.txt", "r");
	//printf("Z jaka dokladnoscia chcesz uzyskac wynik?\n");
	fscanf(fp, "%f", &dokladnosc);
	//printf("%f\n", dokladnosc);

	//printf("\nUstaw punkt startowy ukladu rownan:\n");
	for (i = 0; i < N; i++)
	{
		//printf("x%d = ", i);
		fscanf(fp, "%lf", &x[i]);
		//printf("%lf", x[i]);
	}

	
	printf("Wybierz ile danych chcesz wyswietlac:\n");
	printf("1. Rozwiazania ukladu rownan po kazdym kroku\n");
	printf("2. Samo rozwiazanie rownan\n");
	scanf("%d", &typ);

	//get F(x0)
	getFunction(F, x);

	if (typ != 2)
	{
		/*printf("Wartosc funkcji F: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", F[i]);
		}
		printf("\n");*/
		WriteVectorToFile(F, N, "Wartosc funkcji F: ", counter);
	}




	while ((dnrm2_(&N, dk, &nrhs)>dokladnosc) )// && (counter<=1000)) //norma euklidesowa wektora dk
	{
		
		counter++;

		//F(x)=-F(x)
		dscal_(&N, &alpha, F, &nrhs);

		/*printf("Wartoœci funkcji -F: ");
		for ( i = 0; i < N; i++)
		{
		printf(" %lf ", F[i]);
		}
		printf("\n");*/

		//przypisuje dk=F, zeby dk zostalo nadpisane w mnozeniu i wymnozona zostala wartosc funkcji
		dcopy_(&N, F, &nrhs, dk, &nrhs);

		//1------------------------------------------------------------

		//Bk*dk=-F(x)
		//wyliczam dk
		dgesv_(&N, &nrhs, Bk, &lda, ipiv, dk, &ldb, &info);
		
		//if (info)
			//break;

		/*printf("Macierz dk: ");
		for ( i = 0; i < N; i++)
		{
			printf(" %lf ", dk[i]);
		}
		printf("\n");*/

		//2------------------------------------------------------------------------

		//xk+1=xk+alfa*dk

		//zapisuje starego xk jako dx
		dcopy_(&N, x, &nrhs, dx, &nrhs);

		//xk+1=xk+alfa*dk
		//x staje sie xk+1
		daxpy_(&N, &AlphaK, dk, &nrhs, x, &nrhs);

		/*printf(" xk+1 to: ");
		for ( i = 0; i < N; i++)
		{
			printf(" %lf ", x[i]);
		}
		printf("\n");*/

		//3----------------------------------------------------------------------

		//dx=-dx (gdzie dx to stare xk)
		dscal_(&N, &alpha, dx, &nrhs);

		/*printf("Wektor dx=-dx: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", dx[i]);
		}
		printf("\n");*/

		//y=AlphaK*A*x+beta*y
		//dx=xk+1-xk
		daxpy_(&N, &AlphaK, x, &nrhs, dx, &nrhs);

		/*printf("Wektor dx=xk+1-xk: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", dx[i]);
		}
		printf("\n");*/

		//zapisuje stara funkcje -F(xk)
		dcopy_(&N, F, &nrhs, dF, &nrhs);

		//pobiera funkcje z F(xk+1) x jest teraz xk+1
		getFunction(F, x);
		/*printf("Funkcja xk+1: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", F[i]);
		}
		printf("\n");*/

		//dF=F(xk+1)-F(xk) gdzie -F(xk)=dF
		daxpy_(&N, &AlphaK, F, &nrhs, dF, &nrhs);

		/*printf("dF to: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", dF[i]);
		}
		printf("\n");*/

		//-------------------------------------------------------------------------

		//4 
		//APROKSYMACJA BROYDENA
		//Bk=Bk+(dF-Bdx)*norma(dx)^(-2)*dx**T    dx**T-transponowana macierz dx
		//licznik=dF
		dcopy_(&N, dF, &nrhs, licznik, &nrhs);
		
		//licznik=-Bk*dx+dF
		dgemm_(typeN, typeN, &N, &nrhs, &N, &alpha, Bk, &N, dx, &N, &xAlpha, licznik, &N);
		
		//mianownik=norma(dx)
		//norma wektora
		mianownik[0] = dnrm2_(&N, dx, &nrhs);
		//odwrocenie mianownika w celu mnozenia z licznikiem
		mianownik[0] = 1 / mianownik[0];
		//norma do kwadratu = dx^T*dx
		mianownik[0] = mianownik[0] * mianownik[0];

		//C=alpha*A*B+betaC
		//Bk+1=mianownik*licznik*dx^T+Bk
		dgemm_(typeN, typeT, &N, &N, &nrhs, mianownik, licznik, &N, dx, &N, &xAlpha, Bk, &N);

		/*printf("Macierz Bk:\n");
		for (i = 0; i < N; i++)
		{
			printf("      [");
			for (integer j = i; j < (N*N); j += N)
			{
				printf(" %lf ", Bk[j]);
			}
			printf("]\n");
		}
		printf("\n");*/

	/*	printf("Bk is [%lf %lf %lf]\n", Bk[0], Bk[3], Bk[6]);
		printf("      [%lf %lf %lf]\n", Bk[1], Bk[4], Bk[7]);
		printf("      [%lf %lf %lf]\n", Bk[2], Bk[5], Bk[8]);*/

		//printf("\n\n");

		if (typ != 2)
		{
			/*printf("Rozwiazanie %d to: ", counter);
			for (i = 0; i < N; i++)
			{
				printf("%lf ", x[i]);
			}
			printf("\n");*/
			WriteVectorToFile(x, N, "Rozwiazanie %d to: ", counter);
		}
		

	} 

	//zapis wynikow do pliku
	if (!info) /* succeed */
	{
		printf("\n\nZakonczono z sukcesem!\n\n");
		printf("Rozwiazanie to: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", x[i]);
		}
		printf("\nLiczba krokow: %d \n", counter);
		WriteSolutionToFile(x, N, counter, 1);
	}
	else
	{
		printf("Rozwiazanie to: ");
		for (i = 0; i < N; i++)
		{
			printf("%lf ", x[i]);
		}
		printf("\n");
		fprintf(stderr, "wystapil blad %d\n", info);
		WriteSolutionToFile(x, N, counter, 0);
	}
	
	

	//-------------------------------------------------------------------------

	

	system("pause"); //zatrzymuje okno i czeka na przycisk
	
	clearData();

    return info;
}
