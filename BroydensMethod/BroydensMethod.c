#include < stdio.h>
#include "f2c.h"
#include "clapack.h"
#include "Equations.h"

int main(void)
{
    /* 3x3 matrix A
     * 76 25 11
     * 27 89 51
     * 18 60 32
     */
   /* doublereal A[9] = {76, 27, 18, 25, 89, 60, 11, 51, 32};
    doublereal b[3] = {10, 7, 43};*/

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

	//Bk*dk=-F(x)
	//wyliczam dk
    dgesv_(&N, &nrhs, Bk, &lda, ipiv, dk, &ldb, &info);

	printf("Macierz dk is %lf %lf\n", dk[0], dk[1]);

	//xk+1=xk+alfa*dk
	/*y=AlphaK*A*x+beta*y
	//A - dk, x - Xx, y - xk
	integer Mn = N;
	integer Nn = N;
	char TRANS = 'N';
	
	doublereal beta = 1;
	doublereal Xx[4] = { 1, 0, 0, 1 };

	

	//xk+1=xk+alfa*dk
	//nadpisuje xk jako xk+1
	//dgemv_(&TRANS, &Mn, &Nn, &AlphaK, dk, &lda, Xx, &nrhs, &beta, x, &nrhs);*/
	
	//zapisuje starego x
	dcopy_(&N, x, &nrhs, dx, &nrhs);

	//xk+1=xk+alfa*dk
	daxpy_(&N, &AlphaK, dk, &nrhs, x, &nrhs);

	printf("The xk+1 is %lf %lf\n", x[0], x[1]);

	//dx=-dx (gdzie dx to stare xk)
	dscal_(&N, &alpha, dx, &nrhs);
	
	printf("The dx is %lf %lf\n", dx[0], dx[1]);
	
	//y=AlphaK*A*x+beta*y
	daxpy_(&N, &AlphaK, x, &nrhs, dx, &nrhs);
	//dgemv_(&TRANS, &Mn, &Nn, &AlphaK, x, &lda, Xx, &nrhs, &beta, dx, &nrhs);
	

	printf("The dx is %lf %lf\n", dx[0], dx[1]);
	
	//zapisuje stara funkcje -F(xk)
	dcopy_(&N, F, &nrhs, dF, &nrhs);

	//pobiera funkcje z F(xk+1) x jest teraz xk+1
	getFunction(F, x);
	printf("The xk+1 function is %lf %lf\n", F[0], F[1]);

	//dF=F(xk+1)-F(xk) gdzie -F(xk)=dF
	daxpy_(&N, &AlphaK, F, &nrhs, dF, &nrhs);

	printf("The dF is %lf %lf\n", dF[0], dF[1]);

	if (!info) /* succeed */
		printf("\n\n WELL DONE MALENKA!\n\n");
   else
	fprintf(stderr, "dgesv_ fails %d\n", info);

	system("pause"); //zatrzymuje okno i czeka na przycisk
	//tyrlohg */
    return info;
}
