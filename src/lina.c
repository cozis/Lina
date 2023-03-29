#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include "lina.h"

/* Function: lina_dot
**
**   Evaluates the dot product C = A * B. The A,B
**   matrices are, respectively, mxn and nxl, which
**   means C is mxl. The resulting C matrix is stored
**   in a memory  region specified by the caller. 
**
** Notes:
**
**   - A,B must be provided as contiguous memory regions
**     represented in row-major order. Also, C is stored
**     that way too.
**
**   - The C pointer CAN'T refer to the same memory region 
**     of either A or B.
**
**   - m,n,l must be greater than 0.
**
**   - This function can never fail.
*/
void lina_dot(double *A, double *B, double *C, int m, int n, int l)
{
    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    // Iteration over A's rows
    for(int i = 0; i < m; i++) {

        // Iteration over B's columns
        for(int k = 0; k < l; k++) {

            double sum = 0;

            // Iteration over the single B column 
            // for executing the product of sum
                    
            for(int j=0; j < n; j++)
                sum += A[i * n + j] * B[j * l + k];

            C[i * l + k] = sum;
        }
    }
}

/* Function: lina_add
**
**   Evaluates the matrix addition C = A + B. The result
**   is stored in a memory region provided by the caller.
**   All matrices involved are mxn. 
**
** Notes:
**
**   - A,B must be provided as contiguous memory regions
**     represented in row-major order. Also, C is stored
**     that way too.
**
**   - The C pointer CAN refer to the same memory region 
**     of either A or B.
**
**   - m,n must be greater than 0.
**
**   - This function can never fail.
*/
void lina_add(double *A, double *B, double *C, int m, int n)
{
    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL && C != NULL);
    
    for(int i = 0; i < m*n; i++)
        C[i] = A[i] + B[i];
}

/* Function: lina_scale
**
**   Evaluate B = k * A, where A,B are matrices mxn 
**   and  k is a scalar. The result is stored in a 
**   memory region provided by the caller. 
**
** Notes:
**   - The B pointer CAN refer to the same memory 
**     region of A.
**
**   - m,n must be greater than 0.
**
**   - This function can never fail.
*/
void lina_scale(double *A, double *B, double k, int m, int n){

    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    for(int i = 0; i < m*n; i += 1)
        B[i] = k * A[i];   
}

/* Function: lina_transpose
**
**   Evaluate the transpose of A and store it in B. 
**   The matrix A is mxn, which means B will be nxm.
**
** Notes:
**   - The B pointer CAN refer to the same memory 
**     region of A.
**
**   - m,n must be greater than 0.
**
**   - This function can never fail.
*/
void lina_transpose(double *A, double *B, int m, int n)
{
    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    if(m == 1 || n == 1) {
        // For a matrix with height or width of 1
        // row-major and column-major order coincide,
        // so the stransposition doesn't change the
        // the memory representation. A simple copy
        // does the job.

            if(A != B) // Does the copy or the branch cost more?
                memcpy(B, A, sizeof(A[0]) * m * n);
    
    } else if(m == n) {

        // Iterate over the upper triangular portion of
        // the matrix and switch each element with the
        // corresponding one in the lower triangular portion.
        // NOTE: We're assuming A,B might be the same matrix.
        //       If A,B are the same matrix, then the diagonal
        //       is copied onto itself. By removing the +1 in
        //       the inner loop, the copying of the diagonal
        //       is avoided.

        for(int i = 0; i < n; i += 1)
            for(int j = 0; j < i+1; j += 1) {
                double temp = A[i*n + j];
                B[i*n + j] = A[j*n + i];
                B[j*n + i] = temp;
            }

    } else {
        // Not only the matrix needs to be transposed
        // assuming the destination matrix is the same
        // as the source matrix, but the memory representation
        // of the matrix needs to switch from row-major
        // to col-major, so it's not as simple as switching
        // value's positions.
        // This algorithm starts from the A[0][1] value and
        // moves it where it needs to go, then gets the value
        // that was at that position and puts that in it's
        // new position. This process is iterated until the
        // starting point A[0][1] is overwritten with the
        // new value. In this process the first and last
        // value of the matrix never move.

        B[0] = A[0];
        B[m*n - 1] = A[m*n - 1];

        double item = A[1];
        int    next = m;

        while(next != 1) {
            double temp = A[next];
            B[next] = item;
            item = temp;
            next = (next % n) * m + (next / n);
        }

        B[1] = item;
    }
}

/* Function: scanValue
**
**   Scans a numeric value (such as 12, 4.5, 2.1442) 
**   from the stream [fp] and stores it in [buffer].
**   If more than [max_length] bytes would be written
**   to the buffer, this function fails. The first
**   character of the sequence is assumed to have 
**   been already read and is provided through the
**   [first] argument.
**
**   If the function fails, 0 is returned and an error
**   description is returned through the [error] pointer.
**   If it succeded, then:
**
**     - The [buffer] contains the whole zero-terminated
**       character sequence of the numeric value.
**
**     - Through the [final] pointer is returned the first
**       character that wasn't part of the digit sequence
**       (which was consumed by the function, so if the
**       caller were to read a character from the stream,
**       it would get the second character after the digit
**       sequence).
**
**     - 1 is returned if the sequence represents an integer
**       and -1 if the sequence represents a float.
**   
** Notes:
**   - The buffer is always zero terminated if the
**     function succeded.
**
**   - The [error] and [final] pointers are optional
**     (they can be NULL).
*/
static int scanValue(FILE *fp, char *buffer, int max_length, char first, char *final, char **error)
{
    assert(fp != NULL && buffer != NULL && error != NULL);
    assert(max_length >= 0);
    assert(isdigit(first));

    int n = 0;

    char c = first;

    // Scan the integer portion of
    // the numeric value and copy it
    // into the buffer.
    do {
    
        if(n == max_length) {
            *error = "Internal buffer is too small to hold "
                     "the representation of a numeric value";
            return 0;
        }

        buffer[n++] = c;
            
        c = getc(fp);
    
    } while(c != EOF && isdigit(c));

    // Did the integer part end with
    // a dot?
    _Bool dot = (c == '.');

    // Now scan and copy the decimal 
    // part of the numeric value if 
    // a dot was found.
    if(dot) {
        if(n == max_length) {
            // ERROR: Internal buffer is too small to hold
            //        the representation of this item.
            //        (The dot doesn't fit.)
            *error = "Internal buffer is too small to hold "
                     "the representation of a numeric value";
            return 0;
        }

        buffer[n++] = '.';

        c = getc(fp);

        if(!isdigit(c)) {
            // ERROR: Got something other than a 
            //        digit after the dot.
            *error = "Got something other than a digit after the dot.";
            return 0;
        }
            
        do {
            if(n == max_length) {
                // ERROR: Internal buffer is too small 
                //        to hold the representation of
                //        this item.
                *error = "Internal buffer is too small to hold "
                         "the representation of a numeric value";
                return 0;
            }

            buffer[n++] = c;

            c = getc(fp);
        } while(c != EOF && isdigit(c));
    }

    buffer[n] = '\0';
    
    if(final != NULL)
        *final = c;
    
    return dot ? -1 : 1;
}

/* Function: lina_loadMatrixFromStream
**
**   Load from the stream [fp] a matrix encoded as an
**   ASCII sequence in the form:
**
**     [a b c .. , d e f .. , ..]
**
**   where a,b,c,.. are either integers or floats.
**   For instance, the 4x4 identity matrix is 
**   represented as:
**
**     [1 0 0 0,
**      0 1 0 0,
**      0 0 1 0,
**      0 0 0 1]
**
**   or, equivalently:
**
**     [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]
**
**   since whitespace doesn't matter.
**   The decoded matrix is returned through the return
**   value and is dynamically allocated, therefore the
**   caller must call [free] on it when he doesn't need
**   it anymore. The dimensions of the matrix are returned
**   through the [width] and [height] output arguments.
**
**   If an error occurres (either because an allocation
**   failed or because the matrix syntax is invalid),
**   NULL is returned and a human-readable description of
**   what happened is returned through the [error] pointer.
**
** Notes:
**   - This function skips any whitespace that comes before
**     the matrix in the stream.
**
**   - It can be called multiple times on a stream to get
**     more than one matrix from it.
**
**   - The [error] pointer is optional (it can be NULL).
**
**   - If the stream [fp] is NULL, then [stdin] is used.
*/
double *lina_loadMatrixFromStream(FILE *fp, int *width, int *height, char **error)
{
    assert(width != NULL && height != NULL);
    
    if(fp == NULL)
        fp = stdin;

    char *dummy;
    if(error == NULL)
        error = &dummy;
    else
        *error = NULL;

    char c = getc(fp);

    while(c != EOF && isspace(c))
        c = getc(fp);

    if(c == EOF) {
        // ERROR: Stream ended before a matrix was
        //        found.
        *error = "Stream ended before a matrix was found";
        return NULL;
    }

    if(c != '[') {
        // ERROR: Was expected a '[' as the first 
        //        character of a matrix, but got 
        //        something else instead.
        *error = "Got something other than a matrix "
                 "where one was expected";
        return NULL;
    }

    c = getc(fp);

    // Skip spaces before the first element.
    while(c != EOF && isspace(c))
        c = getc(fp);

    if(c == EOF) {
        // ERROR: Stream ended where a numeric value
        //        was expected. 
        *error = "Stream ended where a numeric value "
                 "was expected";
        return NULL;
    }

    double *matrix = malloc(sizeof(matrix[0]) * 64);

    if(matrix == NULL) {
        // ERROR: Insufficient memory.
        *error = "Insufficient memory";
        return NULL;
    }
    
    int capacity = 64, size = 0, 
        w = -1, i = 0, j = 0;

    if(c != ']')
        while(1) {
            if(!isdigit(c)) {
                // ERROR: Got something other than a digit 
                //        where a numeric value was expected.
                *error = "Got something other than a numeric "
                         "value where one was expected";
                return NULL;
            }

            // Numeric values can't be represented
            // in strings bigger than this buffer
            // since they need to be copied in it
            // to be converted to actual numeric
            // variables.
            char buffer[128];
                
            int res = scanValue(fp, buffer, sizeof(buffer), c, &c, error);

            if(res == 0)
                // Failed to scan the value, abort.
                // NOTE: The error was already reported.
                return NULL;

            assert(res == 1 || res == -1);

            // Make sure the matrix has enough space.
            if(size == capacity) {
                int new_capacity = capacity * 2;

                double *temp = realloc(matrix, sizeof(double) * new_capacity);

                if(temp == NULL) {
                    // ERROR: Insufficient memory.
                    *error = "Insufficient memory";
                    free(matrix);
                    return NULL;
                }

                matrix = temp;
                capacity = new_capacity;
            }

            errno = 0;

            double casted;

            if(res == 1)
                casted = (double) strtoll(buffer, NULL, 10);
            else
                casted = strtod(buffer, NULL);

            if(errno) {
                // ERROR: Failed to convert a numeric value
                //        from it's string form to a numeric 
                //        variable.
                *error = "Failed to convert string to number";
                free(matrix);
                return NULL;
            }

            matrix[size++] = casted;

            i += 1;
                
            while(c != EOF && isspace(c))
                c = getc(fp);

            if(c == ']' || c == ',') {
                // The matrix's row just ended.

                if(w == -1)
                    // This was the first row.
                    w = i;
                else {
                    // This wasn't the first row,
                    // so it's possible that it's
                    // length is different from the
                    // previous ones.
                    assert(w > -1);

                    if(i != w) {
                        // ERROR: The j-th row has the wrong
                        //        number of elements.
                        if(i < w)
                            *error = "Matrix row is too short";
                        else
                            *error = "Matrix row is too long";
                        return NULL;
                    }
                }

                i = 0;
                j += 1;

                if(c == ']')
                    // The whole matrix ended!
                    break;
                    
                c = getc(fp);

                while(c != EOF && isspace(c))
                    c = getc(fp);
            }

            if(c == EOF) {
                // ERROR: Stream ended inside a matrix, where 
                //        either ',', ']' or a numeric value was
                //        expected.
                *error = "Stream ended inside a matrix, where either "
                         "',', ']' or a numeric value was expected";
                return NULL;
            }
        }

    if(size == 0) {
        free(matrix);
        *error = "Empty matrix";
        return NULL;
    }

    // If the internal fragmentation is too much,
    // return a dynamic memory region with the
    // exact size instead of the buffer used to
    // build the matrix.
    int fragm_threshold = 30; // (It's a percentage)
    
    if(100.0 * size/capacity < fragm_threshold) {

        int new_capacity = (size == 0) ? 1 : size;

        double *temp = realloc(matrix, new_capacity * sizeof(double));

        if(temp != NULL)
            matrix = temp;
    }

    *width = w;
    *height = j;

    return matrix;
}

/* Function: lina_saveMatrixToStream
**
**   Save to the stream [fp] a matrix [A] encoding it as an
**   ASCII sequence in the form:
**
**     [a b c .. , d e f .. , ..]
**
**   For instance, the 4x4 identity matrix will
**   be encoded as:
**
**      [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]
**
**   Since the matrix is in row-major order, the caller must
**   specify the collumns and the rows of the matrix
**   through [width] and [height] input arguments.
**
**   If an error occurres, a negative integer is returned 
**   and a human-readable description of what happened 
**   is returned through the [error] pointer.
**
** Notes:
**   - It can be called multiple times on a stream to write
**     more than one matrix on it.
**
**   - The [error] pointer is optional (it can be NULL).
**
**   - If the stream [fp] is NULL, then [stdout] is used.
*/
int lina_saveMatrixToStream(FILE *fp, double *A, int width, int height, char **error)
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
            fprintf(fp, "%f ", A[i*width + j]);

        fprintf(fp, "%f, ", A[i*width + width-1]);
    }

    for (int j = 0; j < width-1; j++)
        fprintf(fp, "%f ", A[(height-1)*width + j]);

    fprintf(fp, "%f", A[(height-1)*width + width-1]);
        
    putc(']',fp);

    return 0;
}

void lina_conv(double *A, double *B, double *C, 
               int Aw, int Ah, int Bw, int Bh)
{
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != B && B != C && C != A);
    assert(Aw > 0 && Ah > 0 && Bw > 0 && Bh > 0);
    assert((Bw & 1) && (Bh & 1)); // B must have odd height and width.

    // NOTE: The output C matrix is smaller than 
    //       A proportionally to B's size.
    
    int Cw = Aw - Bw + 1;
    int Ch = Ah - Bh + 1;
    assert(Cw > 0 && Ch > 0);

    // Iterate over each pixel of the result matrix..
    for(int j = 0; j < Ch; j += 1)
        for(int i = 0; i < Cw; i += 1) {
            // ..and calculate it's value as
            // the scalar product between the
            // mask B and a portion of A.

            C[j * Cw + i] = 0;
            for(int v = 0; v < Bh; v += 1)
                for(int u = 0; u < Bw; u += 1)
                    C[j * Cw + i] += A[(i - Bw/2 + u) * Aw + (i - Bh/2 + v)] * B[v * Bw + u];
        }
}

void lina_reallyP(int *P, double *P2, int n)
{
    memset(P2, 0, sizeof(double) * n * n);

    for (int i = 0; i < n; i++)
        P2[i * n + P[i]] = 1;
}

int lina_decompLUP(double *A, double *L, 
                   double *U, int    *P, 
                   int n)
{
    assert(n > 0);
    assert(A != L && A != U && L != U);

    for (int i = 0; i < n; i++)
        P[i] = i;

    int swaps = 0;
    for (int i = 0; i < n; i++) {

        int v = P[i];
        double max_v = A[v * n + i];
        int    max_i = i;
        
        for (int j = i+1; j < n; j++) {
            int u = P[j];
            double abs = fabs(A[u * n + j]);
            if (abs > max_v) {
                max_v = abs;
                max_i = j;
            }
        }

        if (max_i != i) {

            // Swap rows
            int temp = P[i];
            P[i] = P[max_i];
            P[max_i] = temp;

            swaps++;
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U[i * n + j] = A[P[i] * n + j];

    memset(L, 0, sizeof(double) * n * n);
    for (int i = 0; i < n; i++)
        L[i * n + i] = 1;

    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            double u = U[i * n + i];
            L[j * n + i] = U[j * n + i] / u;
            for (int k = 0; k < n; k++)
                U[j * n + k] -= L[j * n + i] * U[i * n + k];
        }

    return swaps;
}

static void 
printSquareMatrix(double *M, int n, FILE *stream)
{
    for (int i = 0; i < n; i++)
        {
            fprintf(stream, "| ");
            for (int j = 0; j < n; j++)
                {
                    fprintf(stderr, "%2.2f ", M[i * n + j]);
                }
            fprintf(stream, "|\n");
        }
    fprintf(stream, "\n");
}

/* Function: lina_det
**
**   Calculates the determinant of the n by n matrix A
**   and returns it throught the output parameter [det].
**
**   If not enough memory is available, false is returned,
**   else true is returned.
**
** Notes:
**   - The output parameter [det] is optional. (you can
**     ignore the result by passing NULL).
*/
bool lina_det(double *A, int n, double *det)
{
    // Allocate the space for the L,U matrices.
    // I can't think of a version of this algorithm
    // where a temporary buffer isn't necessary.
    double *T = malloc(sizeof(double) * n * n * 2 + sizeof(int) * n);
    if (T == NULL)
        return false;

    // Do the decomposition
    double *L = T;
    double *U = L + (n * n);
    int    *P = (int*) (U + (n * n));
    
    int swaps = lina_decompLUP(A, L, U, P, n);
    if (swaps < 0) {
        free(T);
        return false;
    }

    // Knowing that
    //
    //   A = LU
    //
    // then
    //
    //   det(A) = det(LU) = det(L)det(U)
    //
    // Since L and U are triangular, their 
    // determinant is the product of their 
    // diagonals, so the product of the 
    // determinants is the product of both 
    // the diagonals.

    double prod = 1;
    for (int i = 0; i < n; i++) {
        double l = L[i * n + i];
        double u = U[i * n + i];
        prod *= l * u;
    }

    if (swaps & 1)
        prod = -prod;

    if (det)
        *det = prod;

    free(T);
    return true;
}

/* Checks that [A] is kind of upper triangular.
**
*/
static bool isUpperTriangularEnough(double *A, int n, double eps)
{
    assert(A != NULL && n > 0 && eps > 0);

    // Check that the lower triangular portion (without
    // considering the diagonal) is zero.
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i-1; j++)
            if (A[i * n + j] > eps)
                return false;

    // Now check that the subdiagonal is also zero, 
    // though since we are using the real version of 
    // the QR algorithm, only real eigenvalues can be
    // found. Any comples eigenvalues will manifest 
    // as 2x2 blocks on the diagonal, so we need to 
    // allow such blocks. To do this, a non-zero block 
    // is allowed if it's not following another non-zero
    // block.
    // An important thing to note is that 2x2 matrices
    // will always be considered upper triangular by this
    // function, so the caller must manage this case.
    bool flag = false;
    for (int i = 0; i < n-1; i++) {
        if (fabs(A[(i + 1) * n + i]) > eps) {
            if (flag)
                return false;
            flag = true;
        } else
            flag = false;
    }

    // NOTE: Ideas were taken from [https://math.stackexchange.com/questions/4352389/exact-stop-condition-for-qr-algorithm]
    return true;
}

/* Function: lina_eig
**
**   Calculates the eigenvalues of the n by n matrix M
**   using the (unshifted) QR algorithm and stores them 
**   in the E vector.
**
**   If not enough memory is available, this function
**   aborts returning false. If all went well, true is
**   returned.
**
** Algorithm:
**
**   The algorithm works by decomposing the M matrix into
**   the product of two matrices Q and R, such that Q is
**   orthonormal and R is upper triangular:
**
**       M = QR
**
**   Q and R are then multiplied in inverse order to obtain
**   a new matrix M1, which is then decomposed in two new
**   matrices Q1,R1. The algorithm is iterated n times until 
**   the matrix Mn is upper triangular:
**
**       M = QR  ->  RQ = M(1)
**
**       M(1) = Q(1)R(1)  ->  R(1)Q(1) = M(2)
**
**       M(2) = Q(2)R(2)  ->  R(2)Q(2) = M(3)
**
**              ...
**
**       M(n-1) = Q(n-1)R(n-1)  ->  R(n-1)Q(n-1) = M(n)
**
**       M(n) <--- Triangular!
**
**   The eigenvalues of M(n) are the same as M. Being upper 
**   triangular, M(n) has its eigenvalues on its diagonal,
**   so we just need to scan the diagonal and store it into
**   the E vector. If the original matrix has complex roots,
**   the M(n) sequence will converge to a matrix with a 
**   non-zero 2x2 block on the diagonal for each pair of
**   complex roots. If that's the case, these blocks must
**   be unpacked into the complex values using the quadratic
**   formula.
**   
*/
bool lina_eig(double *M, double complex *E, int n)
{
    // Allocate space for three matrices n by n
    double *T = malloc(sizeof(double) * n * n * 3);
    if (T == NULL)
        return false;

    double *A = T;
    double *Q = A + n * n;
    double *R = Q + n * n;
    memcpy(A, M, sizeof(double) * n * n);

    // At least 100 iterations are done. This is because
    // the QR algorithm doesn't allow complex eigenvalues,
    // so the A matrix may converge to a matrix with 2x2
    // blocks on the diagonal. In general, the algorithm
    // must iterate until the end result is triangular,
    // but that may never be the case, so we end when the
    // result matrix is "kind of triangular" (triangular
    // with 2x2 blocks on the diagonal). But by using this
    // rule, a 2x2 matrix will be considered as tringular
    // from the start, which is not right! That's why we
    // do at least 100 warm-up iterations.
    double eps = 0.1;
    int batch = 100;
    do {
        for (int i = 0; i < batch; i++) {
            lina_decompQR(A, Q, R, n); // A(n) = QR
            lina_dot(R, Q, A, n, n, n); // A(n+1) = RQ
        }
    } while (!isUpperTriangularEnough(A, n, eps));

    // Now we export the diagonal of the iteration result
    // also looking out for 2x2 diagonal blocks, in which
    // case we need to unpack their complex eigenvalues
    for (int i = 0; i < n; i++) {

        // The current diagonal entry is A[i*n + i],
        // so if this is the first entry of a 2x2 block,
        // its lower entry A[(i+1)*n + i] will be non-zero
        if (i+1 < n && fabs(A[(i+1) * n + i]) > eps) {

            // It's a 2x2 block. Unpack the complex eigenvalues
            // using the quadratic formula. (Is there a better
            // way?)
            
            double a = A[(i+0) * n + (i+0)];
            double b = A[(i+0) * n + (i+1)];
            double c = A[(i+1) * n + (i+0)];
            double d = A[(i+1) * n + (i+1)];
            
            // Given the block is:
            //
            //   | a  b |
            //   | c  d |
            //
            // Then the eigenvalues are the roots of:
            //
            //   det(| a-y  b   |) = (a-y)(d-y) - bc = y^2 - (a + d)y + (ad - bc)
            //       | c    d-y |
            //
            // For simplicity:
            //
            //   D = (a + d)^2 - 4(ad - bc)
            //
            // so that
            //
            //   y1, y2 = (a + d)/2 +/- 1/2 sqrt{D}
            //
            // y1 and y2 are one the conjugate of the other. Their
            // real part is
            // 
            //   Re{y1, y2} = (a+d)/2
            //
            // While their immaginary part (in absolute value) is
            //
            //   Imm{y1, y2} = 1/2 sqrt{-D}

            double D = (a+d)*(a+d) - 4*(a*d - b*c);
            assert(D < 0);

            double re = 0.5 * (a+d);
            double im = 0.5 * sqrt(-D);

            double complex y1 = re + im * I;
            double complex y2 = re - im * I;

            // Now place the results into the output vector
            // and tell the loop to skip one iteration
            E[i]   = y1; 
            E[i+1] = y2;
            i++;

        } else
            E[i] = A[i * n + i];
    }

    free(T);
    return true;
}

/* Create the n-1 by n-1 matrix D obtained by
** removing the [del_col] column and [del_row]
** frow the n by n matrix M.
*/
static void 
copyMatrixWithoutRowAndCol(double *M, double *D, int n, 
                           int del_col, int del_row)
{
    // Copy the upper-left portion of matrix M
    // that comes before the deleted column and
    // row.
    for (int i = 0; i < del_row; i++)
        for (int j = 0; j < del_col; j++)
            D[i * (n-1) + j] = M[i * n + j];

    // Copy the lower left portion that comes
    // after both the deleted column and row.
    for (int i = del_row+1; i < n; i++)
        for (int j = del_col+1; j < n; j++)
            D[(i-1) * (n-1) + (j-1)] = M[i * n + j];

    // Copy the bottom portion that comes after
    // the deleted row but before the deleted column.
    for (int i = del_row+1; i < n; i++)
        for (int j = 0; j < del_col; j++)
            D[(i-1) * (n-1) + j] = M[i * n + j];

    // Copy the right portion that comes after
    // the deleted column but before the deleted row.
    for (int i = 0; i < del_row; i++)
        for (int j = del_col+1; j < n; j++)
            D[i * (n-1) + (j-1)] = M[i * n + j];
}

bool lina_inverse(double *M, double *D, int n)
{
    double det;
    if (!lina_det(M, n, &det))
        return false;

    if (det == 0)
        return false; // The matrix can't be inverted

    double *T = malloc(sizeof(double) * ((n-1) * (n-1) + n * n));
    if (T == NULL)
        return false;

    double *M_t = T + (n-1) * (n-1);    
    lina_transpose(M, M_t, n, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            
            copyMatrixWithoutRowAndCol(M_t, T, n, j, i);

            double det2;
            if (!lina_det(T, n-1, &det2)) {
                free(T);
                return false;
            }

            // If the determinant of M isn't zero,
            // neither is this!
            assert(det2 != 0);

            bool i_is_odd = i & 1;
            bool j_is_odd = j & 1;
            int sign = (i_is_odd == j_is_odd) ? 1 : -1;

            D[i * n + j] = sign * det2 / det;
        }

    free(T);
    return true;
}