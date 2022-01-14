
// Evaluate the dot product C = A*B where C is of size m by l, A is m by n and B is n by l.
// NOTE: C can't be the same pointer of A or B.
void lina_dot(double *A, double *B, double *C, int m, int n, int l);

// Evaluate C = A + B where A,B,C are m by n.
// NOTE: C can be the same location of A or B.
void lina_add(double *A, double *B, double *C, int m, int n);

// Evaluate B = k*A where A and B are m by n matrices.
// NOTE: B can be the same location of A.
void lina_scale(double *A, double *B, double k, int m, int n);

// Evaluate B = A^t (the transpose of A) where A is m by n.
// NOTE: B can be the same location of A.
void lina_transpose(double *A, double *B, int m, int n);

// Load from a stream a matrix in the form [a b c.. , d e f.. , ...] 
// where a,b,c,.. are either integer or floating point values.
// The returned pointer must be deallocated using `free`.
double *lina_loadMatrixFromStream(FILE *fp, int *width, int *height, char **error);