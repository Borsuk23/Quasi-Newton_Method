#include "f2c.h"
#include "clapack.h"

/*integer N;
doublereal x[3];
doublereal F[3];
doublereal Bk[9];
doublereal dk[3];
doublereal dx[3];
doublereal dF[3];
integer ipiv[3];*/

integer N;
doublereal* x;
doublereal* F;
doublereal* Bk;
doublereal* Bkpom;
doublereal* dxpom;
doublereal* dk;
doublereal* dx;
doublereal* dF;
doublereal* eye;
integer* ipiv;
doublereal* licznik;
doublereal mianownik[1];
doublereal* work;
float dokladnosc;

void initData();
void initBroyden(doublereal*, doublereal*, float);
void clearData();
getFunction(doublereal*, doublereal*);