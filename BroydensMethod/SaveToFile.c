#include <stdio.h>
#include "SaveToFile.h"



void WriteSolutionsToFile(doublereal* x)
{
	FILE * fp;
	fp = fopen("result.txt", "w");
	
	
	if ((fp = fopen("result.txt", "w")) == NULL)
	{
		printf("Nie mogê otworzyæ pliku do zapisu!\n");
		exit(1);
	}
	else
	{
		fprintf(fp, "The solution is %lf %lf %lf\n", x[0], x[1], x[2]);
		fclose(fp);
	}

	printf("Solutions written to file result.txt\n");
}