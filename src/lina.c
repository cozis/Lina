#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "lina.h"
#include <immintrin.h>
#include <stdint.h>

static void
dot_kernel_6x8(double *A_sub, double *B_sub, double *C_sub, int x, int y, int c_min, int c_max, int n, int l);

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
            // for executing the sum of product
                    
            for(int j=0; j < n; j++)
                sum += A[i * n + j] * B[j * l + k];

            C[i * l + k] = sum;
        }
    }
}

/* Function: lina_dot1
**
**   Evaluates the dot product C = A * B. The A,B
**   matrices are, respectively, mxn and nxl, which
**   means C is mxl. The resulting C matrix is stored
**   in a memory  region specified by the caller. 
**
** Variant 1 of lina_dot:
**   The idea of this variant is that inverting the order
**   of the first and the third loop cicle we can avoid the
**   rolling sum and so breaking the depencency chain 
**   among subsequent add thus increasing the IPC.
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
void lina_dot1(double *A, double *B, double *C, int m, int n, int l)
{
    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    // Since the C matrix can contain any value,
    // this first pass is done to overwrite the values

    // Iteration over A's rows
    for(int i = 0; i < m; i++) {
        // Iteration over B's columns
        for(int k = 0; k < l; k++)
            C[i * l + k] = A[i * n] * B[k];
    }

    // Iteration over the single B column 
    // for executing the sum of product
    for(int j=1; j < n; j++)
    {
        // Iteration over A's rows
        for(int i = 0; i < m; i++) {
            // Iteration over B's columns
            for(int k = 0; k < l; k++)
                C[i * l + k] += A[i * n + j] * B[j * l + k];
        }
    }
}

/* Function: lina_dot2
**
**   Evaluates the dot product C = A * B. The A,B
**   matrices are, respectively, mxn and nxl, which
**   means C is mxl. The resulting C matrix is stored
**   in a memory  region specified by the caller. 
**
** Variant 2 of lina_dot:
**   Other than inverting the order of the first and the
**   third loop cicle this version does the dot product in block
**   of 32x32 values. Doing so the number of cache misses decreases.
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
void lina_dot2(double *A, double *B, double *C, int m, int n, int l)
{
    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    // This size is based on experimental results
    #define BLOCKSIZE 32

    const int br_max = (m & ~(BLOCKSIZE - 1));
    const int bc_max = (l & ~(BLOCKSIZE - 1));

    // Dealing with the squared submatrix of C
    for (int br = 0; br < br_max; br += BLOCKSIZE)
    {
        for (int bc = 0; bc < bc_max; bc += BLOCKSIZE)
        {
            double block[BLOCKSIZE*BLOCKSIZE];
            
            // 1. Compute block

            // Iteration over A's rows
            for(int i = br; i < br+BLOCKSIZE; i++) {

                // Iteration over B's columns
                for(int k = bc; k < bc+BLOCKSIZE; k++)
                    block[(i-br)*BLOCKSIZE + (k-bc)] = A[i * n] * B[k];
            }

            // Iteration over the single B column 
            // for executing the sum of product
            for(int j=1; j < n; j++)
            {
                // Iteration over A's rows
                for(int i = br; i < br+BLOCKSIZE; i++) {

                    // Iteration over B's columns
                    for(int k = bc; k < bc+BLOCKSIZE; k++)
                        block[(i-br)*BLOCKSIZE + (k-bc)] += A[i * n + j] * B[j * l + k];
                }
            }

            // 2. Copy block to C
            for (int i = 0; i < BLOCKSIZE; i++)
                memcpy(&C[(i+br)*l + bc],&block[i*BLOCKSIZE], sizeof(double)*BLOCKSIZE);
        }
    }

    // Dealing with the last rows and cols
    
    // Last rows
    // Iteration over A's rows
    for(int i = br_max; i < m; i++) {
        // Iteration over B's columns
        for(int k = 0; k < l; k++)
            C[i*l + k] = A[i * n ] * B[k];
    }

    // Last cols
    // Iteration over A's rows
    for (int i = 0; i < br_max; i++)
    {
        // Iteration over B's columns
        for(int k = bc_max; k < l; k++)
            C[i*l + k] = A[i * n] * B[k];
    }

    // Iteration over the single B column 
    // for executing the product of sum
    for(int j=1; j < n; j++)
    {
        // Iteration over A's rows
        for(int i = br_max; i < m; i++) {
            // Iteration over B's columns
            for(int k = 0; k < l; k++)
                C[i*l + k] += A[i * n + j] * B[j * l + k];
        }

        // Iteration over A's rows
        for (int i = 0; i < br_max; i++)
        {
            // Iteration over B's columns
            for(int k = bc_max; k < l; k++)
                C[i*l + k] += A[i * n + j] * B[j * l + k];
        }
    }
}

/* Function: lina_dot3
**
**   Evaluates the dot product C = A * B. The A,B
**   matrices are, respectively, mxn and nxl, which
**   means C is mxl. The resulting C matrix is stored
**   in a memory  region specified by the caller. 
**
** Variant 3 of lina_dot:
**   This include the changes of lina_dot2 but uses
**   simd instructions to compute products and sums.
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
void lina_dot3(double *A, double *B, double *C, int m, int n, int l)
{
    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    // This size is based on experimental results
    #define BLOCK_ROWS 6
    #define BLOCK_COLS 8

    const int br_max = (m & ~(BLOCK_ROWS - 1));
    const int bc_max = (l & ~(BLOCK_COLS - 1));

    __m256d *Bm = (__m256d *)B;
    __m256d *Cm = (__m256d *)C;

    // problema: B non è allineato a 32 byte, cosa che pare essere il problema

    // Dealing with the squared submatrix of C
    for (int br = 0; br < br_max; br += BLOCK_ROWS)
    {
        for (int bc = 0; bc < bc_max; bc += BLOCK_COLS)
        {
            __m256d mblock[BLOCK_ROWS][BLOCK_COLS/4] = {0};
            
            // 1. Compute block

            for(int j=0; j < n; j++)
            {
                for(int i = 0; i < BLOCK_ROWS; i++)
                {
                    __m256d A_brdcst = _mm256_broadcast_sd(&A[(i+br) * n + j]);
                    for(int k = 0; k < BLOCK_COLS/4; k++)
                    {
                            mblock[i][k] = _mm256_fmadd_pd(A_brdcst, Bm[(j * l + bc)/4 + k], mblock[i][k]);
                    }
                }
                // Iteration over A's rows
            }

            // 2. Copy block to C
            for (int i = 0; i < BLOCK_ROWS; i++)
                for (int j = 0; j < BLOCK_COLS/4; j++)
                   Cm[((i+br)*l + bc)/4 + j] = mblock[i][j];
        }
    }

    // Dealing with the last rows and cols
    //printf("br_max: %d\nbc_max: %d\n",br_max,bc_max);
    // Last rows
    // Iteration over A's rows
    for(int i = br_max; i < m; i++) {
        // Iteration over B's columns
        for(int k = 0; k < l; k++)
            C[i*l + k] = A[i * n ] * B[k];
    }

    // Last cols
    // Iteration over A's rows
    for (int i = 0; i < br_max; i++)
    {
        // Iteration over B's columns
        for(int k = bc_max; k < l; k++)
            C[i*l + k] = A[i * n] * B[k];
    }

    // Iteration over the single B column 
    // for executing the product of sum
    for(int j=1; j < n; j++)
    {
        // Iteration over A's rows
        for(int i = br_max; i < m; i++) {
            // Iteration over B's columns
            for(int k = 0; k < l; k++)
                C[i*l + k] += A[i * n + j] * B[j * l + k];
        }

        // Iteration over A's rows
        for (int i = 0; i < br_max; i++)
        {
            // Iteration over B's columns
            for(int k = bc_max; k < l; k++)
                C[i*l + k] += A[i * n + j] * B[j * l + k];
        }
    }
}

/* Function: lina_dot4
**
**   Evaluates the dot product C = A * B. The A,B
**   matrices are, respectively, mxn and nxl, which
**   means C is mxl. The resulting C matrix is stored
**   in a memory  region specified by the caller. 
**
** Variant 4 of lina_dot:
**   This include the changes of lina_dot3 but uses the
**   micro kernel subroutine
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
void lina_dot4(double *A, double *B, double *C, int m, int n, int l)
{
    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);
    // A_sub, B_sub and C_sub must be 32 byte aligned
    assert(!((uintptr_t)A & 31llu) && !((uintptr_t)B & 31llu) && !((uintptr_t)C & 31llu));

    #define KERNEL_ROW 6
    #define KERNEL_COLS 8

    const int br_max = (m & ~(KERNEL_ROW - 1));
    const int bc_max = (l & ~(KERNEL_COLS - 1));

    for (int br = 0; br < br_max; br += KERNEL_ROW)
    {
        for (int bc = 0; bc < bc_max; bc += KERNEL_COLS)
        {
            dot_kernel_6x8(A, B, C, br, bc, 0, n, n, l);
        }
    }
}

/*
*   
*   Computes C_sub += A_sub * B_sub where:
*       - C_sub = C[x:x+6][y:y+8]
*       - A_sub = A[x:x+6][c_min:c_max]
*       - B_sub = B[c_min:c_max][y:y+8]
*       - n is the number of columns of A
*       - l the number of columns of B
*/
static void
dot_kernel_6x8(double *A_sub, double *B_sub, double *C_sub, int x, int y, int c_min, int c_max, int n, int l)
{
    // A_sub, B_sub and C_sub must be 32 byte aligned
    // assert is done in the main lina_dot function
    //assert(!((uintptr_t)A_sub & 31llu) && !((uintptr_t)B_sub & 31llu) && !((uintptr_t)C_sub & 31llu));

    // This structure should use 12 YMM registers

    __m256d *Bm_sub = (__m256d *)B_sub;
    __m256d *Cm_sub = (__m256d *)C_sub;
    __m256d acc[6][2] = {0};


    for (int k = c_min; k < c_max; k++)
    {
        for (int i = 0; i < 6; i++)
        {
            __m256d A_brdcst = _mm256_broadcast_sd(&A_sub[(x + i)*n + k]);
            for (int j = 0; j < 2; j++)
                acc[i][j] = _mm256_fmadd_pd(A_brdcst,Bm_sub[(k*l + y)/4 + j],acc[i][j]);
        }
    }

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 2; j++)
            Cm_sub[((x+i)*l + y)/4 + j] = acc[i][j];
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

bool lina_inverse(double *M, double *D, int n)
{
    // To be done
}