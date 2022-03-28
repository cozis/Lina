
/* ---- Operations ---- */
void lina_dot(double *A, double *B, double *C, int m, int n, int l);
void lina_add(double *A, double *B, double *C, int m, int n);
void lina_scale(double *A, double *B, double k, int m, int n);
void lina_transpose(double *A, double *B, int m, int n);
void lina_conv(double *A, double *B, double *C, 
               int Aw, int Ah, int Bw, int Bh);

void lina_dot2(double *A, double *B, double *C, 
               int As, int Bs, int Cs, 
               int m, int n, int l);

/* ---- Utilities ---- */
double *lina_loadMatrixFromStream(FILE *fp, int *width, int *height, char **error);
int     lina_saveMatrixToStream(FILE *fp, double *A, int width, int height, char **error);