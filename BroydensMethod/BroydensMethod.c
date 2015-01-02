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
    doublereal b[3] = {10, 7, 43};

    integer N = 3;
    integer nrhs = 1;
    integer lda = 3;
    integer ipiv[3];
    integer ldb = 3;*/
    integer info;
	doublereal alpha = -1;
	integer nrhs = 1;
	doublereal xAlpha = 1;
	integer lda = 2;


	//get F(x0)
	getFunction(F, x);
	printf("The function is %lf %lf\n", F[0] , F[1] );

	//F(x)=-F(x)
	dscal_(&N, &alpha, F, &nrhs);

	//Bk*dk=-F(x)
   // dgesv_(&N, &nrhs, Bk , &lda, ipiv, b, &ldb, &info);

   if(info) /* succeed */
	printf("The solution is %lf %lf\n", F[0], F[1]);
    else
	fprintf(stderr, "dgesv_ fails %d\n", info);

	system("pause"); //zatrzymuje okno i czeka na przycisk
	//tyrlohg */
    return info;
}
