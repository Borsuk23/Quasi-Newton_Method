#include < stdio.h>
#include "f2c.h"
#include "clapack.h"

int
main(void)
{
    /* 3x3 matrix A
     * 76 25 11
     * 27 89 51
     * 18 60 32
     */
    doublereal A[9] = {76, 27, 18, 25, 89, 60, 11, 51, 32};
    doublereal b[3] = {10, 7, 43};

    integer N = 3;
    integer nrhs = 1;
    integer lda = 3;
    integer ipiv[3];
    integer ldb = 3;
    integer info;
    
    dgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);

    if(info == 0) /* succeed */
	printf("The solution is %lf %lf %lf\n", b[0], b[1], b[2]);
    else
	fprintf(stderr, "dgesv_ fails %d\n", info);

	system("pause"); //zatrzymuje okno i czeka na przycisk

    return info;
}
