#include <stdio.h>
#include "SaveToFile.h"


void WriteSolutionToFile(doublereal* x, int liczbaRownan, int counter, int succes)
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
		if (succes == 1)
		{
			fprintf(fp, "\n\nZakonczono z sukcesem!\n\n");
		}
		fprintf(fp, "Rozwiazanie to: ");
		for (i = 0; i < liczbaRownan; i++)
		{
			fprintf(fp, "%lf ", x[i]);
		}
		fprintf(fp, "\nLiczba krokow: %d \n", counter);

		fclose(fp);

	}
}

void WriteMatrixToFile(doublereal* Bk, int liczbaRownan, int counter)
{
	FILE * fp;
	integer i, j;

	fp = fopen("result.txt", "a");

	if ((fp = fopen("result.txt", "a")) == NULL)
	{
		printf("Nie mogê otworzyæ pliku do zapisu!\n");
		exit(1);
	}
	else
	{
		fprintf(fp, "[%d] Macierz Bk:\n", counter);
		for (i = 0; i < liczbaRownan; i++)
		{
			fprintf(fp, "[");
			for (j = i; j < (liczbaRownan*liczbaRownan); j += liczbaRownan)
			{
				fprintf(fp, "%lf", Bk[j]);
			}
			fprintf(fp, "]\n");
		}
		fprintf(fp, "\n");

		fclose(fp);

	}
}

void WriteVectorToFile(doublereal* f, int liczbaRownan, char* text, int counter)
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
		fprintf(fp, "[%d]", counter);
		fprintf(fp, text);
		for (i = 0; i < liczbaRownan; i++)
		{
			fprintf(fp, "%lf ", f[i]);
		}
		fprintf(fp, "\n");

		fclose(fp);

	}
}

/*void WriteSolutionsToFile(doublereal* x, integer* N)
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

		for (i = 0; i < N; i++)
		{
			fprintf(fp, "x%d = %lf ", i, x[i]);
		}

		fprintf(fp, "\n");

		fclose(fp);
		
	}

	printf("Solutions written to file result.txt\n");
}*/