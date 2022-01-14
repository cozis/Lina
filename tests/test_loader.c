#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "lina.h"

#define check assert

// This program tests the matrix-loading 
// function lina_loadMatrixFromStream.


struct {
	char *src;
	double *matrix;
	int w, h;
} load_tests[] = {
	{ 
		.src = "", 
		.matrix = NULL 
	},
	{ 
		.src = "  []  ", 
		.matrix = NULL 
	},
	{ 
		.src = "  [1]  ", 
		.matrix = (double[]) {1}, 
		.w = 1, 
		.h = 1 
	},
	{ 
		.src = "[1 2 3]", 
		.matrix = (double[]) {1, 2, 3}, 
		.w = 3, 
		.h = 1 
	},
	{ 
		.src = "[1, 2, 3]", 
		.matrix = (double[]) {1, 2, 3}, 
		.w = 1, 
		.h = 3 
	},
	{ 
		.src = "[1 2 3, 4 5 6, 7 8 9]", 
		.matrix = (double[]) {1, 2, 3, 4, 5, 6, 7, 8, 9}, 
		.w = 3, 
		.h = 3 
	},
	{ 
		.src = "[1.0 2.0 3.0, 4.7 5.0 6.0, 7.0 8.0 9.5]", 
		.matrix = (double[]) {1, 2, 3, 4.7, 5, 6, 7, 8, 9.5}, 
		.w = 3, 
		.h = 3 
	},
};

int main()
{
	int total = sizeof(load_tests) / sizeof(load_tests[0]);
	int passed = 0;

	for(int i = 0; i < total; i += 1)
		{
			char *src = load_tests[i].src;

			FILE *stream = fmemopen(src, strlen(src), "r");
			check(stream != NULL);

			int w, h;
			char *err;
			double *M = lina_loadMatrixFromStream(stream, &w, &h, &err);

			double *expected_M = load_tests[i].matrix;

			if(expected_M == NULL)
				{
					// Was expected a failure.
					if(M == NULL)
						{
							// And we got one.
							fprintf(stderr, "Test %d passed.\n", i);
							passed++;
						}
					else
						// Yet the routine succeded.
						fprintf(stderr, "Test %d failed:\n\tWas expected a failure but "
							            "the routine succeded to parse \"%s\"\n", i, src);
				}
			else
				{
					// Was expected a success.
					if(M == NULL)
						// But we got a failure.
						fprintf(stderr, "Test %d failed:\n\tCouldn't parse \"%s\" (%s)\n", i, src, err);
					else
						{
							// And we got a success! Now we need to check
							// that the parsed matrix is the right one.

							int expected_w = load_tests[i].w;
							int expected_h = load_tests[i].h;

							if(expected_w != w || expected_h != h)

								fprintf(stderr, "Test %d failed:\n\tParsing \"%s\" resulted in a "
									    "matrix %dx%d, where a matrix %dx%d was expected\n", 
									    i, src, h, w, expected_h, expected_w);

							else if(memcmp(M, expected_M, sizeof(M[0])*w*h))

								fprintf(stderr, "Test %d failed:\n\tParsing \"%s\" resulted in the "
										"wrong matrix\n", i, src);

							else
								{
									fprintf(stderr, "Test %d passed.\n", i);
									passed += 1;
								}

							free(M);
						}
				}
			fclose(stream);
		}

	fprintf(stderr, "  %d total, %d passed, %d failed\n", total, passed, total - passed);
	return 0;
}