#include <stdio.h>
#include "SaveToFile.h"



void WriteSolutionsToFile(doublereal* x)
{
	FILE * fp;
	integer i;
	fp = fopen("result.txt", "a");
	
	
	if ((fp = fopen("result.txt", "a")) == NULL)
	{
		printf("Nie mogê otworzyæ pliku do zapisu!\n");
		exit(1);
	}
	else
	{
		fprintf(fp, "The solution is: ");

		for (i = 0; i <= 2; i++)
		{
			fprintf(fp, "x%d = %lf ", i, x[i]);
		}

		fprintf(fp, "\n");

		fclose(fp);
		
	}

	printf("Solutions written to file result.txt\n");
}