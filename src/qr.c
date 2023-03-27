#include <math.h>
#include <assert.h>

typedef struct {
    double *items;
    int size;
} square_matrix_t;

typedef struct {
    double *items;
    int stride;
    int size;
} vector_t;

typedef struct {
    vector_t base;
    double scale;
} scaled_vector_t;

static square_matrix_t 
square_matrix_from_raw(double *M, int n)
{
    return (square_matrix_t) {.items=M, .size=n};
}

static vector_t 
get_column_of_square_matrix(square_matrix_t M, int i)
{
    assert(i > -1 && i < M.size);

    return (vector_t) { 
        .items  = M.items + i, 
        .stride = M.size, 
        .size   = M.size
    };
}

static void 
copy_vector(vector_t V, vector_t S)
{
    assert(V.size == S.size);
    for (int i = 0; i < V.size; i++)
        V.items[V.stride * i] = S.items[S.stride * i];
}

static void 
subtract_vector_inplace(vector_t V, scaled_vector_t S)
{
    assert(V.size == S.base.size);

    for (int i = 0; i < V.size; i++)
        V.items[V.stride * i] -= S.scale * S.base.items[S.base.stride * i];
}

static void
scale_vector_inplace(vector_t V, double a)
{
    for (int i = 0; i < V.size; i++)
        V.items[V.stride * i] *= a;
}

static scaled_vector_t
scale_vector_lazily(vector_t V, double a)
{
    return (scaled_vector_t) {.base=V, .scale=a};
}

static double 
scalar_product(vector_t V, vector_t U)
{
    assert(V.size == U.size);

    double scale = 0;
    for (int i = 0; i < V.size; i++)
        scale += V.items[i * V.stride] * U.items[i * U.stride];
    return scale;
}

static double 
calculate_norm(vector_t V)
{
    double sum_of_squares = scalar_product(V, V);
    return sqrt(sum_of_squares);
}

static double 
normalize_inplace(vector_t V)
{
    double norm = calculate_norm(V);
    if (norm != 0)
        scale_vector_inplace(V, 1/norm);
    return norm;
}

static scaled_vector_t
project(vector_t V, vector_t U)
{
    double scale_vu = scalar_product(V, U);
    double scale_uu = scalar_product(U, U);
    double ratio = scale_vu / scale_uu;
    return scale_vector_lazily(U, ratio);
}

/** Gram-Schmidt orthonormalization
 **/
void lina_orthoNormGramSchmidt(double *A, double *Q, int n)
{
    square_matrix_t A2 = square_matrix_from_raw(A, n);
    square_matrix_t Q2 = square_matrix_from_raw(Q, n);

    for (int i = 0; i < n; i++) {

        vector_t Qi = get_column_of_square_matrix(Q2, i);
        vector_t Ai = get_column_of_square_matrix(A2, i);
        copy_vector(Qi, Ai);

        for (int j = 0; j < i; j++) {
            vector_t Qj = get_column_of_square_matrix(Q2, j);
            subtract_vector_inplace(Qi, project(Ai, Qj));
        }

        normalize_inplace(Qi);
        // TODO: Handle case of zero norm 
    }
}

void lina_decompQR(double *A, double *Q, double *R, int n)
{
    lina_orthoNormGramSchmidt(A, Q, n);

    square_matrix_t A2 = square_matrix_from_raw(A, n);
    square_matrix_t Q2 = square_matrix_from_raw(Q, n);

    // Now calculate R by multiplying Q^t and A
    for(int i = 0; i < n; i++) { // Iterate over each column i of Q..
        for(int j = 0; j < n; j++) { // ..and over each column j of A
            vector_t Qi = get_column_of_square_matrix(Q2, i);
            vector_t Aj = get_column_of_square_matrix(A2, j);
            R[i * n + j] = scalar_product(Qi, Aj);
        }
    }
}