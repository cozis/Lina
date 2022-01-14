#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "lina.h"

/* This function compares the lina_transpose
** implementation against the naive implementation.
*/

#define check assert

static void naive_transpose(double *A, double *B, int m, int n)
{
    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    double *support = malloc(sizeof(*support) * m * n);

    check(support != NULL);

    memcpy(support, A, sizeof(*support) * m * n);

    for(int i = 0; i < n; i++)
        
        for(int j = 0; j < m; j++)

            B[j*n + i] = support[i*m + j];

    free(support);
}

int main()
{
	int m = 1000;
	int n = 100000;

	double *big = malloc(sizeof(double) * m * n);
	check(big != NULL);

	memset(big, 0, sizeof(double) * m * n);

	double t1, t2;
	clock_t begin, end;

	begin = clock();
	lina_transpose(big, big, m, n);
	end = clock();
	t1 = (double) (end - begin) / CLOCKS_PER_SEC;

	begin = clock();
	naive_transpose(big, big, m, n);
	end = clock();
	t2 = (double) (end - begin) / CLOCKS_PER_SEC;

	printf("lina_transpose took %gms\n", t1*1000);
	printf("naive_transpose took %gms\n", t2*1000);

	free(big);
	return 0;
}