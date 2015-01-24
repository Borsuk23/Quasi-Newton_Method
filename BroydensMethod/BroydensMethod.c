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
	integer lwork = 2 * N;
	integer sizeB = N*N;
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

	

	//pierwsza macierz Broydena
	initBroyden(F, x, dokladnosc);

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




	while ((dnrm2_(&N, dk, &nrhs)>dokladnosc)  && (counter<=10000000)) //norma euklidesowa wektora dk
	{
		
		counter++;

		//F(x)=-F(x)
		//dscal_(&N, &alpha, F, &nrhs);

		//zapisuje F do dk
		dcopy_(&N, F, &nrhs, dk, &nrhs);
		//zapisuje Bk do Bkpom
		dcopy_(&sizeB, Bk, &nrhs, Bkpom, &nrhs);
		//skaluje dk*-1
		dscal_(&N, &alpha, dk, &nrhs);

		

		//1------------------------------------------------------------
		//Rozwiazanie Bkdk=-F dk sie nadpisuje, bo wczesniej bylo F
		//Bk*dk=-F(x)
		//wyliczam dk
		dgesv_(&N, &nrhs, Bk, &N, ipiv, dk, &N, &info);

		
		if (info)
			break;

		/*printf("Macierz dk: ");
		for ( i = 0; i < N; i++)
		{
			printf(" %lf ", dk[i]);
		}
		printf("\n");*/

		//2------------------------------------------------------------------------

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
		//Obliczanie dF=Fk+1
		getFunction(dF, x);

		//zapisanie do pamieci zeby potem przekazac do nastepnej iteracji
		//zapisuje dF do dx
		dcopy_(&N, dF, &nrhs, dx, &nrhs);

		// obliczenie dF staje sie y y=-F+dF
		daxpy_(&N, &alpha, F, &nrhs, dF, &nrhs);

		

		//-------------------------------------------------------------------------

		//4 
		//APROKSYMACJA BROYDENA
		//Bk=Bk+(dF-Bdx)*norma(dx)^(-2)*dx**T    dx**T-transponowana macierz dx
		//y=dF y staje sie licznikiem
		dcopy_(&N, dF, &nrhs, licznik, &nrhs);
		//wiemy ze dk to roznica xow wiec zamiast dx jest dk
		//-Bkpom*dx+y staje sie licznikiem
		dgemm_(typeN, typeN, &N, &nrhs, &N, &alpha, Bkpom, &N, dk, &N, &xAlpha, licznik, &N);
		//licznik*dx^T+zero*Bk staje sie Bk
		dgemm_(typeN, typeT, &N, &N, &nrhs, &AlphaK, licznik, &N, dk, &N, &zero, Bk, &N);
		//dx^T*dx+zero*mianownik staje sie mianownikiem
		dgemm_(typeT, typeN, &nrhs, &nrhs, &N, &AlphaK, dk, &N, dk, &N, &zero, mianownik, &nrhs);
		mianownik[0] = 1 / mianownik[0];
		//mianownik*Bk(czyli licznik)*dxpom + Bkpom (poprzednia macierz Broydena) staje sie nowym Bkpom
		dgemm_(typeN, typeN, &N, &N, &N, mianownik, Bk, &N, dxpom, &N, &AlphaK, Bkpom, &N);
		//zapisz Bkpom do Bk
		dcopy_(&sizeB, Bkpom, &nrhs, Bk, &nrhs);
		//zapisz dx do F (dzieki temu jest przepisane Fk+1 do nowego Fk)
		dcopy_(&N, dx, &nrhs, F, &nrhs);

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
		WriteMatrixToFile(Bk, N, counter);
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
		WriteMatrixToFile(Bk, N, counter);
	}
	
	

	//-------------------------------------------------------------------------

	

	system("pause"); //zatrzymuje okno i czeka na przycisk
	
	clearData();

    return info;
}
