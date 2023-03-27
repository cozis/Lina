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

void print_vector(double *V, int n, FILE *stream)
{
    fprintf(stream, "[ ");
    for (int i = 0; i < n; i++)
        fprintf(stderr, "%2.2f ", V[i]);
    fprintf(stream, "]\n");
}

int main(void)
{
    double A[4] = {1, 2, 3, 4};
    
    double L[4];
    double U[4];
    double LU[4];
    lina_decompLU(A, L, U, 2);
    lina_dot(L, U, LU, 2, 2, 2);
    print_square_matrix(A,  2, stderr);
    print_square_matrix(L,  2, stderr);
    print_square_matrix(U,  2, stderr);
    print_square_matrix(LU, 2, stderr);

    double det;
    lina_det(A, 2, &det);
    fprintf(stderr, "det=%2.2f\n", det);

    double Q[4];
    double R[4];
    double QR[4];
    lina_decompQR(A, Q, R, 2);
    lina_dot(Q, R, QR, 2, 2, 2);
    print_square_matrix(Q,  2, stderr);
    print_square_matrix(R,  2, stderr);
    print_square_matrix(QR, 2, stderr);

    double E[4];
    lina_eig(A, E, 2);
    print_vector(E, 2, stderr);

    return 0;
}
