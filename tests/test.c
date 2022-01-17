#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "lina.h"

#define check assert

//Print the matrix A with size m by n
static void pmatrix(FILE *fp, double *A, int m, int n);

struct dot_test{

	double *A;
	double *B;
	double *C;

};



int main()
{
	struct dot_test a;
	

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