#include < stdio.h>
#include "f2c.h"
#include "clapack.h"
#include "Equations.h"
#include "SaveToFile.h"

int main(void)
{
    
    integer info;
	float dokladnosc;
	doublereal alpha = -1;
	doublereal zero = 0;
	integer nrhs = 1;
	doublereal xAlpha = 1;
	integer lda = N;
	integer ldb = N;
	doublereal AlphaK = 1;
	integer counter = 0;
	char typeN[1] = { 'N' }; //macierz normalna
	char typeT[1] = { 'T' }; //macierz transponowana
	

	initData();

	//pobor z konsoli dokladnosci wybranej przez uzytkownika
	printf("Z jaka dokladnoscia chcesz uzyskac wynik?");
	scanf("%f", &dokladnosc);

	//pobor z konsoli punktu startowego wybranego przez uzytkownika
	printf("\nUstaw punkt startowy ukladu rownan:\n");
	for (integer i = 0; i <= N-1 ; i++)
	{
		printf("x%d\n", i);
		scanf("%lf", &x[i]);
	}


	//get F(x0)
	getFunction(F, x);
	printf("The x0 function is \n");
	for (integer i = 0; i <= N - 1; i++)
	{
		printf(" %lf\n", F[i]);
	}
	
	


	while ((dnrm2_(&N, dk, &nrhs)>dokladnosc) && (counter<=1000)) //norma euklidesowa wektora dk
	{
		
		counter++;

		//F(x)=-F(x)
		dscal_(&N, &alpha, F, &nrhs);

		printf("The -F function is\n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", F[i]);
		}

		//przypisuje dk=F, zeby dk zostalo nadpisane w mnozeniu i wymnozona zostala wartosc funkcji
		dcopy_(&N, F, &nrhs, dk, &nrhs);

		//1------------------------------------------------------------

		//Bk*dk=-F(x)
		//wyliczam dk
		dgesv_(&N, &nrhs, Bk, &lda, ipiv, dk, &ldb, &info);
		
		if (info)
			break;

		printf("Macierz dk is \n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", dk[i]);
		}
		//2------------------------------------------------------------------------

		//xk+1=xk+alfa*dk

		//zapisuje starego xk jako dx
		dcopy_(&N, x, &nrhs, dx, &nrhs);

		//xk+1=xk+alfa*dk
		//x staje sie xk+1
		daxpy_(&N, &AlphaK, dk, &nrhs, x, &nrhs);

		printf("The xk+1 is \n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", x[i]);
		}
		//3----------------------------------------------------------------------

		//dx=-dx (gdzie dx to stare xk)
		dscal_(&N, &alpha, dx, &nrhs);

		printf("The dx (-dx) is \n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", dx[i]);
		}

		//y=AlphaK*A*x+beta*y
		//dx=xk+1-xk
		daxpy_(&N, &AlphaK, x, &nrhs, dx, &nrhs);

		printf("The dx is \n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", dx[i]);
		}

		//zapisuje stara funkcje -F(xk)
		dcopy_(&N, F, &nrhs, dF, &nrhs);

		//pobiera funkcje z F(xk+1) x jest teraz xk+1
		getFunction(F, x);
		printf("The xk+1 function is \n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", F[i]);
		}
		//dF=F(xk+1)-F(xk) gdzie -F(xk)=dF
		daxpy_(&N, &AlphaK, F, &nrhs, dF, &nrhs);

		printf("The dF is \n");
		for (integer i = 0; i <= N - 1; i++)
		{
			printf(" %lf\n", dF[i]);
		}
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

		

		printf("Bk is [%lf %lf %lf]\n", Bk[0], Bk[3], Bk[6]);
		printf("      [%lf %lf %lf]\n", Bk[1], Bk[4], Bk[7]);
		printf("      [%lf %lf %lf]\n", Bk[2], Bk[5], Bk[8]);

		printf("\n\n");


	} 
	
	if (!info) /* succeed */
	{
		printf("\n\n WELL DONE !\n\n");
		printf("The solution is %lf %lf %lf\n", x[0], x[1], x[2]);
		printf("The counter is %li\n", counter);
	}
	else
		fprintf(stderr, "dgesv_ fails %d\n", info);

	//-------------------------------------------------------------------------

	//zapis wyników do pliku
	WriteSolutionsToFile(x);

	system("pause"); //zatrzymuje okno i czeka na przycisk
	
	clearData();

    return info;
}
