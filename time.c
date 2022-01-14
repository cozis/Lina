#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "lina.h"

/* This program compares the lina_transpose
** implementation against the naive implementation.
*/

#define check assert

static void naive_transpose(double *A, double *B, int m, int n)
{
    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    double *support;

    if(A == B)
    	{
    		support = malloc(sizeof(*support) * m * n);
    		
    		check(support != NULL);

		    memcpy(support, A, sizeof(*support) * m * n);
    	}
    else
    	{
    		support = A;
    	}

    for(int i = 0; i < n; i++)
        
        for(int j = 0; j < m; j++)

            B[j*n + i] = support[i*m + j];

    if(support != A)
	    free(support);
}

// Wrap transposing functions and return their
// execution time.
static double time_transposition(void (*callback)(double*, double*, int, int), double *A, double *B, int m, int n)
{
	clock_t begin = clock();

	callback(A, B, m, n);

	clock_t end = clock();
	
	return (double) (end - begin) / CLOCKS_PER_SEC;
}

int main()
{
	int m = 1000;
	int n = 100000;

	double *big = malloc(sizeof(double) * m * n);
	check(big != NULL);

	memset(big, 0, sizeof(double) * m * n);

	printf("lina_transpose took %gms (in-place)\n", 
		1000 * time_transposition(lina_transpose, big, big, m, n));

	printf("naive_transpose took %gms (in-place)\n", 
		1000 * time_transposition(naive_transpose, big, big, m, n));

	double *big2 = malloc(sizeof(double) * m * n);
	check(big2 != NULL);

	printf("lina_transpose took %gms\n", 
		1000 * time_transposition(lina_transpose, big, big2, m, n));
	
	printf("naive_transpose took %gms\n", 
		1000 * time_transposition(naive_transpose, big, big2, m, n));

	free(big);
	free(big2);
	return 0;
}