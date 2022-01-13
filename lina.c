#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "lina.h"

#define check assert

void lina_dot(double *A, double *B, double *C, int m, int n, int l){

    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    //Actual dot
    //Iteration over A's rows
    for(int i=0; i < m;i++){

        //Iteration over B's columns
        for(int k=0; k < l; k++){

            double pos = 0;
            
            //Iteration over the single B column for executing the product of sum
            for(int j=0; j < n; j++)
                pos += A[i * n + j] * B[j * l + k];

            //The usage of a support variable 'pos' is for safety purpose (non-zero C matrix)
            C[i*l+k] = pos;

        }
        
    }

}

void lina_add(double *A, double *B, double *C, int m, int n){

    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL && C != NULL);
    
    for(int i = 0; i < m; i++)

        for(int j = 0; j < n; j++)

            C[i * n + j] = A[i * n + j] + B[i * n + j];
}

void lina_scale(double *A, double *B, double k, int m, int n){

    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    for(int i = 0; i < m; i++)

        for(int j = 0; j < n; j++)

            B[i * n + j] = k * A[i * n + j];
}

void lina_transpose(double *A, double *B, int m, int n){

    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

#warning "Try to think better algorithm"

    double *support = malloc(sizeof(*support) * m * n);

    check(support != NULL);

    memcpy(support, A, sizeof(*support) * m * n);

    for(int i=0;i<n;i++){

        for(int j=0;j<m;j++){

            B[i*m + j] = support[i*m + j];
        }
    
    }

    free(support);
}