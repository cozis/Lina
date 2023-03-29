#include <complex.h>
#include <stdbool.h>

void   lina_dot(double *A, double *B, double *C, int m, int n, int l);
void   lina_add(double *A, double *B, double *C, int m, int n);
bool   lina_det(double *A, int n, double *det);
void   lina_scale(double *A, double *B, double k, int m, int n);
void   lina_transpose(double *A, double *B, int m, int n);
bool   lina_inverse(double *M, double *D, int n);
void   lina_conv(double *A, double *B, double *C, 
                 int Aw, int Ah, int Bw, int Bh);

void lina_dot2(double *A, double *B, double *C, 
               int As, int Bs, int Cs, 
               int m, int n, int l);

bool lina_eig(double *M, double complex *E, int n);
void lina_reallyP(int *P, double *P2, int n);
int  lina_decompLUP(double *A, double *L, double *U, int *P, int n);
void lina_decompQR(double *A, double *Q, double *R, int n);
void lina_orthoNormGramSchmidt(double *A, double *Q, int n);

double *lina_loadMatrixFromStream(FILE *fp, int *width, int *height, char **error);
int     lina_saveMatrixToStream(FILE *fp, double *A, int width, int height, char **error);