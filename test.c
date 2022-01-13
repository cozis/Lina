#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "lina.h"

#define NUM_TESTS 100           //Number of tests case to create
#define MATRIX_MAX_ROWS 4       //Max rows that a matrix can have
#define MATRIX_MAX_COLLUMNS 4   //Max collumns a matrix can have
#define MAX_VALUE 15            //Max value a matrix can have
//Matrix type. A is m by n

typedef struct matrix_t
{
	double *A;
	int m;
	int n;
}matrix_t;

//Dinamically create NUM_TESTS tests case randomly generating matrices with random values (whom max is MAX_VALUE)
//and random collumns and row numbers (whom max is MATRIX_MAX_ROWS and MATRIX_MAX_COLLUMNS)
void create_tests(matrix_t **tests);
//Free the heap memory from the test structure
void delete_tests(matrix_t **tests);

//Print the matrix A with size m by n
void pmatrix(double *A, int m, int n);

int main()
{
	
	
	matrix_t *tests[NUM_TESTS];

	create_tests(tests);

	pmatrix(tests[0]->A,tests[0]->m,tests[0]->n);

	return 0;
}


void pmatrix(double *A,int m,int n){

	for(int i = 0; i<m; i++){

		for(int j = 0; j< n; j++)
			printf("%f\t",A[i*n + j]);
		
		printf("\n");

	}


}

void create_tests(matrix_t **tests){

	srand(1);

	for(int i=0;i<NUM_TESTS;i++){

		int rows = rand() % MATRIX_MAX_ROWS;
		int cols = rand() % MATRIX_MAX_COLLUMNS;
		
		double *ptr = malloc(sizeof(*ptr)*rows*cols);
		
		for(int j=0;j<rows*cols;j++){

			ptr[j] = rand() % MAX_VALUE;

		}

		matrix_t *m_point = malloc(sizeof(matrix_t));

		m_point->A = ptr;
		m_point->m = rows;
		m_point->n = cols;

		tests[i] = m_point;
	}
	

}

void delete_tests(matrix_t **tests){

	for(int i=0;i<NUM_TESTS;i++){

		free(tests[i]->A);
		free(tests[i]);

	}

	tests = NULL;

}
