
// Evaluate the dot product C = A*B where C is of size m by l, A is m by n and B is n by l.
// NOTE: C can't be the same pointer of A or B.
void lina_dot(double *A, double *B, double *C, unsigned int m, unsigned int n, unsigned int l);

// Evaluate C = A + B where A,B,C are m by n.
// NOTE: C can be the same location of A or B.
void lina_add(double *A, double *B, double *C, unsigned int m, unsigned int n);

// Evaluate B = k*A where A and B are m by n matrices.
// NOTE: B can be the same location of A.
void lina_scale(double *A, double *B, double k, unsigned int m, unsigned int n);

// Evaluate B = A^t (the transpose of A) where A is m by n.
// NOTE: B can be the same location of A.
void lina_transpose(double *A, double *B, unsigned int m, unsigned int n);