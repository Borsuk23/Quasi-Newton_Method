#include < stdio.h>
#include "f2c.h"
#include "clapack.h"
#include "Equations.h"

int main(void)
{
    
    integer info;
	doublereal alpha = -1;
	doublereal zero = 0;
	integer nrhs = 1;
	doublereal xAlpha = 1;
	integer lda = N;
	integer ldb = N;
	doublereal dx[3] = { 0, 0, 0 };
	doublereal AlphaK = 1;
	//integer counter = 0;
	char typeN[1] = { 'N' };
	char typeT[1] = { 'T' };
	doublereal licznik[3] = { 0, 0, 0 };
	doublereal mianownik[1] = { 1 };
	doublereal mianownik_pom[3] = { 1, 1, 1 };
	doublereal eye[1] = { 1 };
	doublereal mianownik_odwrotnosc[1] = { 1 };
	doublereal broyden_pom[3] = { 1, 1, 1 };

	
	//get F(x0)
	getFunction(F, x);
	printf("The function is %lf %lf %lf\n", F[0], F[1], F[2]);

	while (dnrm2_(&N, dk, &nrhs)>0.0001) //norma euklidesowa wektora dk
	{
		

		//F(x)=-F(x)
		dscal_(&N, &alpha, F, &nrhs);

		printf("The -F function is %lf %lf %lf\n", F[0], F[1], F[2]);

		//przypisuje dk=F, zeby dk zostalo nadpisane w mnozeniu i wymnozona zostala wartosc funkcji
		dcopy_(&N, F, &nrhs, dk, &nrhs);

		//1------------------------------------------------------------

		//Bk*dk=-F(x)
		//wyliczam dk
		dgesv_(&N, &nrhs, Bk, &lda, ipiv, dk, &ldb, &info);

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
		//Bk liczone jako Bk^(-1)
		//licznik=dx
		dcopy_(&N, dx, &nrhs, licznik, &nrhs);
		//licznik=-Bk*dF+dx
		dgemm_(typeN, typeN, &N, &nrhs, &N, &alpha, Bk, &N, dF, &N, &xAlpha, licznik, &N);
		//mianownik_pom=dx**T*B
		dgemm_(typeT, typeN, &nrhs, &N, &N, &xAlpha, dx, &N, Bk, &N, &zero, mianownik_pom, &nrhs);
		//mianownik=mianownik_pom*mianownik
		dgemm_(typeN, typeN, &nrhs, &nrhs, &N, &xAlpha, mianownik_pom, &nrhs, dF, &N, &zero, mianownik, &nrhs);
		//wyznaczenie mianownik^-1 za pomoca Ax=B, gdzie B=eye, wtedy x=A^(-1)
		//mianownik_odwrotnosc=eye
		dcopy_(&nrhs, eye, &nrhs, mianownik_odwrotnosc, &nrhs);
		//mianownik_odwrotnosc=mianownik^(-1)
		dgesv_(&nrhs, &nrhs, mianownik, &nrhs, ipiv, mianownik_odwrotnosc, &nrhs, &info);
		//broyden_pom=licznik*mianownik_odwrotnosc
		dgemm_(typeN, typeN, &N, &nrhs, &nrhs, &xAlpha, licznik, &N, mianownik_odwrotnosc, &nrhs, &zero, broyden_pom, &N);
		//Bk=Bk+broyden_pom*mianownik_pom
		dgemm_(typeN, typeN, &N, &N, &nrhs, &xAlpha, broyden_pom, &N, mianownik_pom, &nrhs, &xAlpha, Bk, &N);

		printf("Bk is [%lf %lf %lf]\n", Bk[0], Bk[3], Bk[6]);
		printf("      [%lf %lf %lf]\n", Bk[1], Bk[4], Bk[7]);
		printf("      [%lf %lf %lf]\n", Bk[2], Bk[5], Bk[8]);


		if (!info) /* succeed */
		{
			printf("\n\n WELL DONE MALENKA !\n\n");
			//printf("counter %li\n\n", counter);
		}
		else
			fprintf(stderr, "dgesv_ fails %d\n", info);

	} 
	
	
	//-------------------------------------------------------------------------

	

	system("pause"); //zatrzymuje okno i czeka na przycisk
	
    return info;
}
