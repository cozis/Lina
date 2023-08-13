#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/lina.h"

#define A_ROWS 1000llu
#define A_COLS 146llu

#define B_ROWS 146llu
#define B_COLS 1024llu

uint64_t nanos();

int main()
{
    uint64_t ops = A_ROWS*B_COLS*2*A_COLS;
    uint64_t start,stop,lina_dot_time, lina_dot_mod1_time, lina_dot_mod1_1_time, lina_dot_mod2_time;
    double *A = (double *)malloc(sizeof(double)*A_ROWS*A_COLS);
    double *B = (double *)malloc(sizeof(double)*B_ROWS*B_COLS);

    double *C1 = (double *)malloc(sizeof(double)*A_ROWS*B_COLS);
    double *C2 = (double *)malloc(sizeof(double)*A_ROWS*B_COLS);
    double *C3 = (double *)malloc(sizeof(double)*A_ROWS*B_COLS);
    double *C4 = (double *)malloc(sizeof(double)*A_ROWS*B_COLS);
    
    for (int i = 0; i < A_ROWS*A_COLS; i++)
        A[i] = (double)(rand()%10);
    for (int i = 0; i < B_ROWS*B_COLS; i++)
        B[i] = (double)(rand()%10);
    for (int i = 0; i < A_ROWS*B_COLS; i++)
    {
        C1[i] = (double)(rand()%10);
        C2[i] = (double)(rand()%10);
        C3[i] = (double)(rand()%10);
        C4[i] = (double)(rand()%10);
    }

    start = nanos();
    lina_dot(A,B,C1,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot_time = stop-start;
    
    start = nanos();
    lina_dot_mod1(A,B,C2,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot_mod1_time = stop-start;

    start = nanos();
    lina_dot_mod2(A,B,C3,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot_mod2_time = stop-start;

    start = nanos();
    lina_dot_mod2_old(A,B,C4,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot_mod1_1_time = stop-start;

    if(!memcmp(C1,C2,sizeof(double)*A_ROWS*B_COLS) && !memcmp(C2,C3,sizeof(double)*A_ROWS*B_COLS) && !memcmp(C3,C4,sizeof(double)*A_ROWS*B_COLS))
    {
        printf( "lina_dot     : %f GFLOPS\n"
                "lina_dot_mod1: %f GFLOPS\n"
                "lina_dot_mod2: %f GFLOPS\n"
                "lina_dot_mod2_old: %f GFLOPS\n", (double)ops/lina_dot_time,
                                                (double)ops/lina_dot_mod1_time,
                                                (double)ops/lina_dot_mod2_time,
                                                (double)ops/lina_dot_mod1_1_time);
        
    }
    else
        printf("ERRORE: i prodotti matriciali sono diversi!\n");
    
    free(A);
    free(B);
    free(C1);
    free(C2);
    free(C3);
    free(C4);
}