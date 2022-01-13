#include <stddef.h> // NULL
#include <assert.h> // assert
#include "lina.h"

// Evaluate the dot product C = A*B where C is of size m by l, A is m by n and B is n by l.
// NOTE: C can't be the same pointer of A or B.
void lina_dot(double *A, double *B, double *C, int m, int n, int l){

    assert(m >= 0 && n >= 0 && l >= 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    //Actual dot
    //Iteration on A's rows
    for(int i=0; i < m;i++){
        
        double pos = 0;

        //Iteration on B's collumns
        for(int k=0; k < l; k++){

            //Iteration on the single B collumn for executing the product of sum
            for(int j=0; j < n; j++)
                pos += A[i * n + j] * B[j * l + k];

            //The usage of a support variable 'pos' is for safety purpose (non-zero C matrix)
            C[i*l+k] = pos;

        }
        
    }

}

// Evaluate C = A + B where A,B,C are m by n.
// NOTE: C can be the same location of A or B.
void lina_add(double *A, double *B, double *C, int m, int n){

    assert(m >= 0 && n >= 0);
    assert(A != NULL && B != NULL && C != NULL);
    
    for(int i = 0; i < m; i++)

        for(int j = 0; j < n; j++)

            C[i * n + j] = A[i * n + j] + B[i * n + j];
}

// Evaluate B = k*A where A and B are m by n matrices.
// NOTE: B can be the same location of A.
void lina_scale(double *A, double *B, double k, int m, int n){

    assert(m >= 0 && n >= 0);
    assert(A != NULL && B != NULL);

    for(int i = 0; i < m; i++)

        for(int j = 0; j < n; j++)

            B[i * n + j] = k * A[i * n + j];
}

// Evaluate B = A^t (the transpose of A) where A is m by n.
// NOTE: B can be the same location of A.
void lina_transpose(double *A, double *B, int m, int n){

}