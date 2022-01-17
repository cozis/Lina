#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "lina.h"

#define check assert

//Print the matrix A with size m by n
static void pmatrix(FILE *fp, double *A, int m, int n);



int main()
{	
	int n,m;
	double *M;
	char *msg;
	

	FILE *ptr = fopen("tests/dot/t1.txt","r");
	check(ptr != NULL);	

	M = lina_loadMatrixFromStream(ptr, &n, &m, &msg);

	printf("%d, %d", m,n);
	//pmatrix(stdout,M,m,n);

	free(M);
	fclose(ptr);
	return 0;
}

static void pmatrix(FILE *fp, double *A,int m,int n){

	for(int i = 0; i<m; i++){
		fprintf(fp, "    | ");
		for(int j = 0; j< n; j++)
			fprintf(fp, "%g ",A[i*n + j]);
		
		fprintf(fp, "|\n");

	}

	fprintf(fp, "\n");

}