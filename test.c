#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "lina.h"

#define NUM_TESTS 6      //Number of tests case to create
#define MAX_VALUE 15       //Max value a matrix can have

typedef struct matrix_t
{
	double *A;
	int m;
	int n;
}matrix_t;

//Dinamically create NUM_TESTS tests case randomly generating matrices with random values (whom max is MAX_VALUE)
//index from 0 to NUM_TESTS/3 will contain 2x3 matrices, index from NUM_TESTS/3 to 2*NUM_TESTS/3 3x3, and index from 2*NUM_TESTS/3 to NUM_TESTS 4x3
void create_tests(matrix_t **tests);

//Free the heap memory from the test structure
void delete_tests(matrix_t **tests);

//Print the matrix A with size m by n
void pmatrix(double *A, int m, int n);

int main()
{
	matrix_t *tests[NUM_TESTS];

	create_tests(tests);

	for(int i=0;i< NUM_TESTS;i++)
		pmatrix(tests[i]->A,tests[i]->m,tests[i]->n);

	double *C = malloc(sizeof(*C)*tests[0]->m*tests[0]->n);

	lina_scale(tests[0]->A,C,1.5,tests[0]->m,tests[0]->n);

	pmatrix(C,tests[0]->m,tests[0]->n);

	delete_tests(tests);
	return 0;
}

void pmatrix(double *A,int m,int n){

	for(int i = 0; i<m; i++){

		for(int j = 0; j< n; j++)
			printf("%f\t",A[i*n + j]);
		
		printf("\n");

	}

	printf("\n");

}

void create_tests(matrix_t **tests){

	srand(1);

	for(int i=0;i<NUM_TESTS;i++){

		int rows;
		int cols;

		if(i < NUM_TESTS/3){
			rows = 2;
			cols = 3;
		}
		else if (i >= NUM_TESTS/3 && i < 2*NUM_TESTS/3)
		{
			rows = 3;
			cols = 3;
		}
		else{
			rows = 4;
			cols = 3;
		}
		

		
		
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
