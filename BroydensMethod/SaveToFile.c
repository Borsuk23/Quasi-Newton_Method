#include <stdio.h>
#include "SaveToFile.h"



void WriteSolutionsToFile(doublereal* x)
{
	FILE * fp;
	integer i;
	fp = fopen("result.txt", "w");
	
	
	if ((fp = fopen("result.txt", "w")) == NULL)
	{
		printf("Nie mogê otworzyæ pliku do zapisu!\n");
		exit(1);
	}
	else
	{
		fprintf(fp, "The solution is:\n");

		for (i = 0; i <= 2; i++)
		{
			fprintf(fp, "x%d = %lf\n", i, x[i]);
		}
		fclose(fp);
		
	}

	printf("Solutions written to file result.txt\n");
}