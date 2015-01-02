#include < stdio.h>
#include "f2c.h"
#include "clapack.h"
#include "Equations.h"

int main(void)
{
    
    integer info;
	doublereal alpha = -1;
	integer nrhs = 1;
	doublereal xAlpha = 1;
	integer lda = N;
	integer ldb = N;
	doublereal dx[2] = { 0, 0 };
	doublereal AlphaK = 1;

	//get F(x0)
	getFunction(F, x);
	printf("The function is %lf %lf\n", F[0] , F[1] );

	//F(x)=-F(x)
	dscal_(&N, &alpha, F, &nrhs);

	printf("The -F function is %lf %lf\n", F[0], F[1]);

	//przypisuje dk=F, zeby dk zostalo nadpisane w mnozeniu i wymnozona zostala wartosc funkcji
	dcopy_(&N, F, &nrhs, dk, &nrhs);

	//1------------------------------------------------------------

	//Bk*dk=-F(x)
	//wyliczam dk
    dgesv_(&N, &nrhs, Bk, &lda, ipiv, dk, &ldb, &info);

	printf("Macierz dk is %lf %lf\n", dk[0], dk[1]);

	//2------------------------------------------------------------------------

	//xk+1=xk+alfa*dk
	
	//zapisuje starego xk jako dx
	dcopy_(&N, x, &nrhs, dx, &nrhs);

	//xk+1=xk+alfa*dk
	//x staje sie xk+1
	daxpy_(&N, &AlphaK, dk, &nrhs, x, &nrhs);

	printf("The xk+1 is %lf %lf\n", x[0], x[1]);

	//3----------------------------------------------------------------------

	//dx=-dx (gdzie dx to stare xk)
	dscal_(&N, &alpha, dx, &nrhs);
	
	printf("The dx is %lf %lf\n", dx[0], dx[1]);
	
	//y=AlphaK*A*x+beta*y
	//dx=xk+1-xk
	daxpy_(&N, &AlphaK, x, &nrhs, dx, &nrhs);
	
	printf("The dx is %lf %lf\n", dx[0], dx[1]);
	
	//zapisuje stara funkcje -F(xk)
	dcopy_(&N, F, &nrhs, dF, &nrhs);

	//pobiera funkcje z F(xk+1) x jest teraz xk+1
	getFunction(F, x);
	printf("The xk+1 function is %lf %lf\n", F[0], F[1]);

	//dF=F(xk+1)-F(xk) gdzie -F(xk)=dF
	daxpy_(&N, &AlphaK, F, &nrhs, dF, &nrhs);

	printf("The dF is %lf %lf\n", dF[0], dF[1]);

	//-------------------------------------------------------------------------

	if (!info) /* succeed */
		printf("\n\n WELL DONE MALENKA !\n\n");
   else
	fprintf(stderr, "dgesv_ fails %d\n", info);

	system("pause"); //zatrzymuje okno i czeka na przycisk
	
    return info;
}
