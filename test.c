#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "lina.h"

#define check assert

//Print the matrix A with size m by n
static void pmatrix(FILE *fp, double *A, int m, int n);

struct {
	double A[9], // Left argument
		   B[9], // Right argument
		   C[9]; // Expected result
} dot_tests[] = {
	{ 
		.A = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
		},
		.B = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
		},
		.C = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
		},
	},
};

struct {
	double A[9], // Input
	       B[9], // Expected output
	       factor; // Scale factor
} scale_tests[] = {
	{
		.A = {
			1, 1, 1,
			1, 1, 1,
			1, 1, 1,
		},
		.B = {
			2, 2, 2,
			2, 2, 2,
			2, 2, 2,
		},
		.factor = 2,
	},
};

struct {
	double *A, *B;
	int m, n;
} transp_tests[] = {
	{
		.A = (double[]) {
			1, 0, 0,
			0, 2, 0,
			0, 0, 3,
		},
		.B = (double[]) {
			1, 0, 0,
			0, 2, 0,
			0, 0, 3,
		},
		.m = 3,
		.n = 3,
	},
	{
		.A = (double[]) {
			1, 2, 3,
			0, 0, 0,
			0, 0, 0,
		},
		.B = (double[]) {
			1, 0, 0,
			2, 0, 0,
			3, 0, 0,
		},
		.m = 3,
		.n = 3,
	},
	{
		.A = (double[]) {
			0, 1, 0,
			0, 2, 0,
			0, 3, 0,
		},
		.B = (double[]) {
			0, 0, 0,
			1, 2, 3,
			0, 0, 0,
		},
		.m = 3,
		.n = 3,
	},
	{
		.A = (double[]) {
			0, 0, 1,
			0, 0, 2,
			0, 0, 3,
		},
		.B = (double[]) {
			0, 0, 0,
			0, 0, 0,
			1, 2, 3,
		},
		.m = 3,
		.n = 3,
	},
	{
		.A = (double[]) {
			1, 2, 3,
			4, 5, 6,
		},
		.B = (double[]) {
			1, 4,
			2, 5,
			3, 6,
		},
		.m = 2,
		.n = 3,
	},
	{
		.A = (double[]) {
			1, 4,
			2, 5,
			3, 6,
		},
		.B = (double[]) {
			1, 2, 3,
			4, 5, 6,
		},
		.m = 3,
		.n = 2,
	},
};

int main()
{
	// Evaluate dot product tests.
	{
		int dot_passed = 0;
		int dot_total  = sizeof(dot_tests) / sizeof(*dot_tests);

		for(int i = 0; i < dot_total; i += 1)
			{
				double R[9];

				lina_dot(dot_tests[i].A, dot_tests[i].B, R, 3, 3, 3);

				if(!memcmp(R, dot_tests[i].C, sizeof(R)))
					{
						fprintf(stderr, "Dot product test %d passed.\n", i);
						dot_passed += 1;
					}
				else
					{
						fprintf(stderr, "Dot product test %d failed:\n  got matrix:\n\n", i);
						pmatrix(stderr, R, 3, 3);
						fprintf(stderr, "  instead of:\n\n");
						pmatrix(stderr, dot_tests[i].C, 3, 3);
					}
			}
		fprintf(stderr, "\n\t%d dot products out of %d were succesful.\n\n", dot_passed, dot_total);
	}

	
	// Evaluate scaling tests.
	{
		int scale_passed = 0;
		int scale_total  = sizeof(scale_tests) / sizeof(*scale_tests);

		for(int i = 0; i < scale_total; i += 1)
			{
				double R[9];

				lina_scale(scale_tests[i].A, R, scale_tests[i].factor, 3, 3);

				if(!memcmp(R, scale_tests[i].B, sizeof(R)))
					{
						fprintf(stderr, "Scaling test %d passed.\n", i);
						scale_passed += 1;
					}
				else
					{
						fprintf(stderr, "Scaling test %d failed:\n  got matrix:\n\n", i);
						pmatrix(stderr, R, 3, 3);
						fprintf(stderr, "  instead of:\n\n");
						pmatrix(stderr, scale_tests[i].B, 3, 3);
					}
			}
		fprintf(stderr, "\n\t%d scalings out of %d were succesful.\n\n", scale_passed, scale_total);
	}

	// Evaluate transposition tests.
	{
		int transp_passed = 0;
		int transp_total  = sizeof(transp_tests) / sizeof(*transp_tests);

		for(int i = 0; i < transp_total; i += 1)
			{
				int m = transp_tests[i].m;
				int n = transp_tests[i].n;

				double R[32];

				assert(sizeof(R) >= m * n * sizeof(R[0]));

				lina_transpose(transp_tests[i].A, R, m, n);

				if(!memcmp(R, transp_tests[i].B, sizeof(R[0]) * m * n))
					{
						fprintf(stderr, "Transposition test %d passed.\n", i);
						transp_passed += 1;
					}
				else
					{
						fprintf(stderr, "Transposition test %d failed:\n  got matrix:\n\n", i);
						pmatrix(stderr, R, m, n);
						fprintf(stderr, "  instead of:\n\n");
						pmatrix(stderr, transp_tests[i].B, m, n);
					}
			}
		fprintf(stderr, "\n\t%d transpositions out of %d were succesful.\n\n", transp_passed, transp_total);
	}

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