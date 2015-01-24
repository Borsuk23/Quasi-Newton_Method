#include "f2c.h"
#include "clapack.h"


void WriteSolutionToFile(doublereal* x, int liczbaRownan, int counter, int succes);
void WriteMatrixToFile(doublereal* Bk, int liczbaRownan, int counter);
void WriteVectorToFile(doublereal* f, int liczbaRownan, char* text, int counter);