#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../src/lina.h"

#define A_ROWS 960llu
#define A_COLS 960llu

#define B_ROWS 960llu
#define B_COLS 960llu

int saveMatrixToStream(FILE *fp, double *A, int width, int height, char **error);
static uint64_t nanos();

int main()
{
    uint64_t ops = A_ROWS*B_COLS*2*A_COLS;
    uint64_t start,stop,lina_dot_time, lina_dot1_time, lina_dot2_time, lina_dot3_time, lina_dot4_time;
    
    double *A = (double *)aligned_alloc(32,sizeof(double)*A_ROWS*A_COLS);
    double *B = (double *)aligned_alloc(32,sizeof(double)*B_ROWS*B_COLS);

    double *C1 = (double *)aligned_alloc(32,sizeof(double)*A_ROWS*B_COLS);
    double *C2 = (double *)aligned_alloc(32,sizeof(double)*A_ROWS*B_COLS);
    double *C3 = (double *)aligned_alloc(32,sizeof(double)*A_ROWS*B_COLS);
    double *C4 = (double *)aligned_alloc(32,sizeof(double)*A_ROWS*B_COLS);
    double *C5 = (double *)aligned_alloc(32,sizeof(double)*A_ROWS*B_COLS);
    
    for (int i = 0; i < A_ROWS*A_COLS; i++)
        A[i] = (double)(rand()%2);
    for (int i = 0; i < B_ROWS*B_COLS; i++)
        B[i] = (double)(rand()%2);
    for (int i = 0; i < A_ROWS*B_COLS; i++)
    {
        C1[i] = (double)(rand()%2);
        C2[i] = (double)(rand()%2);
        C3[i] = (double)(rand()%2);
        C4[i] = (double)(rand()%2);
        C5[i] = (double)(rand()%2);
    }

    start = nanos();
    lina_dot(A,B,C1,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot_time = stop-start;
    
    start = nanos();
    lina_dot1(A,B,C2,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot1_time = stop-start;

    start = nanos();
    lina_dot2(A,B,C3,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot2_time = stop-start;

    start = nanos();
    lina_dot3(A,B,C4,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot3_time = stop-start;

    start = nanos();
    lina_dot4(A,B,C5,A_ROWS,A_COLS,B_COLS);
    stop = nanos();

    lina_dot4_time = stop-start;

    if(!memcmp(C1,C2,sizeof(double)*A_ROWS*B_COLS)
        && !memcmp(C2,C3,sizeof(double)*A_ROWS*B_COLS)
        && !memcmp(C3,C4,sizeof(double)*A_ROWS*B_COLS)
        && !memcmp(C4,C5,sizeof(double)*A_ROWS*B_COLS))
    {
        printf( "lina_dot : %f GFLOPS\n"
                "lina_dot1: %f GFLOPS\n"
                "lina_dot2: %f GFLOPS\n"
                "lina_dot3: %f GFLOPS\n"
                "lina_dot4: %f GFLOPS\n",   (double)ops/lina_dot_time,
                                            (double)ops/lina_dot1_time,
                                            (double)ops/lina_dot2_time,
                                            (double)ops/lina_dot3_time,
                                            (double)ops/lina_dot4_time);
        FILE *fp = fopen("lina_dots_success.txt", "w");
        if (!fp)
            return -1;
        
        saveMatrixToStream(fp,C1,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C2,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C3,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C4,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C5,A_ROWS,A_COLS,NULL);

        fclose(fp);
    }
    else
    {
        printf("ERRORE: i prodotti matriciali sono diversi!\n");
        FILE *fp = fopen("lina_dots_error.txt", "w");
        if (!fp)
            return -1;
        
        saveMatrixToStream(fp,C1,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C2,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C3,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C4,A_ROWS,A_COLS,NULL);
        fprintf(fp,"\nFINE\n");
        saveMatrixToStream(fp,C5,A_ROWS,A_COLS,NULL);

        fclose(fp);
    }
    
    free(A);
    free(B);
    free(C1);
    free(C2);
    free(C3);
    free(C4);
    free(C5);
}

static uint64_t nanos()
{
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return (uint64_t)time.tv_sec*1000000000 + (uint64_t)time.tv_nsec;
}

int saveMatrixToStream(FILE *fp, double *A, int width, int height, char **error)
{
    assert(A != NULL);

    char *dummy;
    if (error == NULL)
        error = &dummy;
    else
        *error = NULL;
    
    if (width < 1) {
        *error = "The provided width is less than one";
        return -1;
    }
    
    if (height < 1) {
        *error = "The provided height is less than one";
        return -1;
    }

    if (fp == NULL)
        fp = stdout;

    putc('[',fp);

    for (int i = 0; i < height-1; i++) {
        for (int j = 0; j < width-1; j++)
            fprintf(fp, "%.1f ", A[i*width + j]);

        fprintf(fp, "%.1f,\n", A[i*width + width-1]);
    }

    for (int j = 0; j < width-1; j++)
        fprintf(fp, "%.1f ", A[(height-1)*width + j]);

    fprintf(fp, "%.1f", A[(height-1)*width + width-1]);
        
    putc(']',fp);

    return 0;
}