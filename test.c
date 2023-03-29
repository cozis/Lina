#include <stdio.h>
#include "src/lina.h"

void print_square_matrix(double *M, int n, FILE *stream)
{
    for (int i = 0; i < n; i++)
        {
            fprintf(stream, "| ");
            for (int j = 0; j < n; j++)
                {
                    fprintf(stderr, "%2.2f ", M[i * n + j]);
                }
            fprintf(stream, "|\n");
        }
    fprintf(stream, "\n");
}

void print_vector(double complex *V, int n, FILE *stream)
{
    fprintf(stream, "[ ");
    for (int i = 0; i < n; i++)
        fprintf(stderr, "(%2.2f + i%2.2f) ", creal(V[i]), cimag(V[i]));
    fprintf(stream, "]\n");
}

int main(void)
{
    
    double M[25] = {
        1, 2, 3, 4, 5,
        5, 1, 2, 3, 4,
        4, 5, 1, 2, 3,
        3, 4, 5, 1, 2,
        2, 3, 4, 5, 1,
    };

    fprintf(stderr, "# --- M --- #\n");
    print_square_matrix(M,  5, stderr);
    

    /*
    double  L[25];
    double  U[25];
    int     P[5];
    double P2[25];
    lina_decompLUP(M, L, U, P, 5);
    lina_reallyP(P, P2, 5);

    fprintf(stderr, "# --- L --- #\n");
    print_square_matrix(L, 5, stderr);

    fprintf(stderr, "# --- U --- #\n");
    print_square_matrix(U, 5, stderr);

    fprintf(stderr, "# --- P2 --- #\n");
    print_square_matrix(P2, 5, stderr);

    double PA[25];
    lina_dot(P2, M, PA, 5, 5, 5);
    fprintf(stderr, "# --- PA --- #\n");
    print_square_matrix(PA,  5, stderr);

    double LU[25];
    lina_dot(L, U, LU, 5, 5, 5);
    fprintf(stderr, "# --- LU --- #\n");
    print_square_matrix(LU, 5, stderr);

    double det;
    lina_det(M, 5, &det);
    fprintf(stderr, "det=%2.2f\n", det);

    fprintf(stderr, "# --- eig(M) --- #\n");
    double complex E[5];
    lina_eig(M, E, 5);
    print_vector(E, 5, stderr);
    */

    
    double invM[25];
    lina_inverse(M, invM, 5);

    double expI[25];
    lina_dot(M, invM, expI, 5, 5, 5);

    fprintf(stderr, "# --- inv(M) --- #\n");
    print_square_matrix(invM, 5, stderr);

    fprintf(stderr, "# --- I? --- #\n");
    print_square_matrix(expI, 5, stderr);
    

    /*
    double M[16] = {
        1, 5, 4, 2,
        2, 1, 5, 3,
        4, 3, 2, 5,
        5, 4, 3, 1,
    };

    fprintf(stderr, "# --- M --- #\n");
    print_square_matrix(M, 4, stderr);

    double L[16];
    double U[16];
    int P[4];
    lina_decompLUP(M, L, U, P, 4);
    fprintf(stderr, "# --- L,U,P --- #\n");
    print_square_matrix(L, 4, stderr);
    print_square_matrix(U, 4, stderr);

    fprintf(stderr, "[ ");
    for(int i = 0; i < 4; i++)
        fprintf(stderr, "%d ", P[i]);
    fprintf(stderr, "]\n");

    double RP[16];
    double PM[16];
    lina_reallyP(P, RP, 4);
    lina_dot(RP, M, PM, 4, 4, 4);

    double LU[16];
    lina_dot(L, U, LU, 4, 4, 4);

    fprintf(stderr, "# --- P,PM,LU --- #\n");
    print_square_matrix(RP, 4, stderr);
    print_square_matrix(PM, 4, stderr);
    print_square_matrix(LU, 4, stderr);

    double det;
    lina_det(M, 4, &det);

    fprintf(stderr, "det(M) = %2.2f\n", det);
    */
    return 0;
}
