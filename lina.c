#include "lina.h"

// Evaluate the dot product C = A*B where C is of size m by l, A is m by n and B is n by l.
// NOTE: C can't be the same pointer of A or B.
void lina_dot(double **A, double **B, double **C, unsigned int m, unsigned int n, unsigned int l){

    //Control for the sanity of the inputs:
    //C must point to a different memory than A and B
    if(C == A || C == B)
        //Implement some Cozis-like-good-way method for error signaling
        return;

    //Actual dot
    //Iteration on A's rows
    for(int i=0; i < m;i++){
        
        double pos = 0;

        //Iteration on B's collumns
        for(int k=0; k < l; k++){

            //Iteration on the single B collumn for executing the product of sum
            for(int j=0; j < n; j++){

            pos = pos + A[i][j] * B[j][k];

            }

            //The usage of a support variable 'pos' is for safety purpose (non-zero C matrix)
            C[i][k] = pos;

        }
        
    }

}

// Evaluate C = A + B where A,B,C are m by n.
// NOTE: C can be the same location of A or B.
void lina_add(double **A, double **B, double **C, unsigned int m, unsigned int n){

    for(int i=0;i<m;i++){

        for(int j=0;j<n;j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }

}

// Evaluate B = k*A where A and B are m by n matrices.
// NOTE: B can be the same location of A.
void lina_scale(double **A, double **B, double k, unsigned int m, unsigned int n){

    for(int i=0;i<m;i++){

        for(int j=0;j<n;j++){
            B[i][j] = k * A[i][j];
        }
    }

}

// Evaluate B = A^t (the transpose of A) where A is m by n.
// NOTE: B can be the same location of A.
void lina_transpose(double **A, double **B, unsigned int m, unsigned int n){

    

}