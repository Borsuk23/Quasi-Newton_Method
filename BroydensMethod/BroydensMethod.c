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
	doublereal dx[3] = { 0, 0, 0 };
	doublereal AlphaK = 1;
	integer counter = 0;
	char typeN[1] = { 'N' };
	char typeT[1] = { 'T' };
	doublereal licznik[3] = { 0, 0, 0 };
	doublereal mianownik[1] = { 1 };
	doublereal mianownik_pom[3] = { 1, 1, 1 };
	doublereal eye[1] = { 1 };
	doublereal mianownik_odwrotnosc[1] = { 1 };
	doublereal broyden_pom[3] = { 1, 1, 1 };

	//pobor z konsoli dokladnosci wybranej przez uzytkownika
	printf("Z jaka dokladnoscia chcesz uzyskac wynik?");
	scanf("%f", &dokladnosc);

	printf("\nUstaw punkt startowy ukladu rownan:\n");
	for (integer i = 0; i <= 2; i++)
	{
		printf("x%d\n", i);
		scanf("%lf", &x[i]);
	}

	//get F(x0)
	getFunction(F, x);
	printf("The function is %lf %lf %lf\n", F[0], F[1], F[2]);

	while ((dnrm2_(&N, dk, &nrhs)>dokladnosc) && (counter<=1000)) //norma euklidesowa wektora dk
	{
		
		counter++;

		//F(x)=-F(x)
		dscal_(&N, &alpha, F, &nrhs);

		printf("The -F function is %lf %lf %lf\n", F[0], F[1], F[2]);

		//przypisuje dk=F, zeby dk zostalo nadpisane w mnozeniu i wymnozona zostala wartosc funkcji
		dcopy_(&N, F, &nrhs, dk, &nrhs);

		//1------------------------------------------------------------

		//Bk*dk=-F(x)
		//wyliczam dk
		dgesv_(&N, &nrhs, Bk, &lda, ipiv, dk, &ldb, &info);
		
		if (info)
			break;

		printf("Macierz dk is %lf %lf %lf\n", dk[0], dk[1], dk[2]);

		//2------------------------------------------------------------------------

		//xk+1=xk+alfa*dk

		//zapisuje starego xk jako dx
		dcopy_(&N, x, &nrhs, dx, &nrhs);

		//xk+1=xk+alfa*dk
		//x staje sie xk+1
		daxpy_(&N, &AlphaK, dk, &nrhs, x, &nrhs);

		printf("The xk+1 is %lf %lf %lf\n", x[0], x[1], x[2]);

		//3----------------------------------------------------------------------

		//dx=-dx (gdzie dx to stare xk)
		dscal_(&N, &alpha, dx, &nrhs);

		printf("The dx is %lf %lf %lf\n", dx[0], dx[1], dx[2]);

		//y=AlphaK*A*x+beta*y
		//dx=xk+1-xk
		daxpy_(&N, &AlphaK, x, &nrhs, dx, &nrhs);

		printf("The dx is %lf %lf %lf\n", dx[0], dx[1], dx[2]);

		//zapisuje stara funkcje -F(xk)
		dcopy_(&N, F, &nrhs, dF, &nrhs);

		//pobiera funkcje z F(xk+1) x jest teraz xk+1
		getFunction(F, x);
		printf("The xk+1 function is %lf %lf %lf\n", F[0], F[1], F[2]);

		//dF=F(xk+1)-F(xk) gdzie -F(xk)=dF
		daxpy_(&N, &AlphaK, F, &nrhs, dF, &nrhs);

		printf("The dF is %lf %lf %lf\n", dF[0], dF[1], dF[2]);

		//-------------------------------------------------------------------------

		//4 
		//Aproksymacja Broydena
		//licznik=dF
		dcopy_(&N, dF, &nrhs, licznik, &nrhs);
		
		//licznik=-Bk*dx+dF
		dgemm_(typeN, typeN, &N, &nrhs, &N, &alpha, Bk, &N, dx, &N, &xAlpha, licznik, &N);
		
		//norma wektora
		mianownik[0] = dnrm2_(&N, dx, &nrhs);
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
	
    return info;
}
