#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include "lina.h"

#define check assert

void lina_dot(double *A, double *B, double *C, int m, int n, int l){

    assert(m > 0 && n > 0 && l > 0);
    assert(A != NULL && B != NULL && C != NULL);
    assert(A != C && B != C);

    // Iteration over A's rows
    for(int i = 0; i < m; i++)
        {
            // Iteration over B's columns
            for(int k = 0; k < l; k++)
                {
                    double pos = 0;

                    // Iteration over the single B column 
                    // for executing the product of sum
                    
                    for(int j=0; j < n; j++)

                        pos += A[i*n + j] * B[j*l + k];

                    C[i*l + k] = pos;
                }
        }
}

void lina_add(double *A, double *B, double *C, int m, int n){

    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL && C != NULL);
    
    for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)

            C[i*n + j] = A[i*n + j] + B[i*n + j];
}

void lina_scale(double *A, double *B, double k, int m, int n){

    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    for(int i = 0; i < m; i++)

        for(int j = 0; j < n; j++)

            B[i*n + j] = k * A[i*n + j];
}

void lina_transpose(double *A, double *B, int m, int n)
{
    assert(m > 0 && n > 0);
    assert(A != NULL && B != NULL);

    if(m == 1 || n == 1)
        {
            memcpy(B, A, sizeof(A[0]) * m * n);
            return;
        }

    if(m == n)
        {
            for(int i = 0; i < n; i += 1)
                for(int j = 0; j < i+1; j += 1)
                    {
                        double temp = A[i*n + j];
                        B[i*n + j] = A[j*n + i];
                        B[j*n + i] = temp;
                    }
        }
    else
        {
            B[0] = A[0];
            B[m*n-1] = A[m*n-1];

            double item = A[1];
            int    next = m;

            while(next != 1)
                {
                    double temp = A[next];
                    B[next] = item;
                    item = temp;
                    next = (next % n) * m + (next / n);
                }

            B[1] = item;
        }
}

// Returns 0 if an error occurred, 1 if an integer 
// was scanned and -1 if a floating point was scanned.
static int scanValue(FILE *fp, char *buffer, int max_length, char first, char *final, char **error)
{
    assert(fp != NULL && buffer != NULL && error != NULL);
    assert(isdigit(first));

    int n = 0;

    char c = first;

    // Scan the integer portion of
    // the numeric value and copy it
    // into the buffer.
    do
        {
            if(n == max_length)
                {
                    // ERROR: Internal buffer is too small to hold
                    //        the representation of this item.
                    *error = "Internal buffer is too small to hold "
                             "the representation of a numeric value";
                    return 0;
                }

            buffer[n++] = c;
            
            c = getc(fp);
        }
    while(c != EOF && isdigit(c));

    // Did the integer part end with
    // a dot?
    _Bool dot = (c == '.');

    // Now scan and copy the decimal 
    // part of the numeric value if 
    // a dot was found.
    if(dot)
        {
            if(n == max_length)
                {
                    // ERROR: Internal buffer is too small to hold
                    //        the representation of this item.
                    //        (The dot doesn't fit.)
                    *error = "Internal buffer is too small to hold "
                             "the representation of a numeric value";
                    return 0;
                }

            buffer[n++] = '.';

            c = getc(fp);

            if(!isdigit(c))
                {
                    // ERROR: Got something other than a 
                    //        digit after the dot.
                    *error = "Got something other than a digit after the dot.";
                    return 0;
                }
            
            do
                {
                    if(n == max_length)
                        {
                            // ERROR: Internal buffer is too small 
                            //        to hold the representation of
                            //        this item.
                            *error = "Internal buffer is too small to hold "
                                     "the representation of a numeric value";
                            return 0;
                        }

                    buffer[n++] = c;

                    c = getc(fp);
                }
            while(c != EOF && isdigit(c));
        }

    buffer[n] = '\0';
    
    if(final != NULL)
        *final = c;
    
    return dot ? -1 : 1;
}

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

    if(c == EOF)
        {
            // ERROR: Stream ended before a matrix was
            //        found.
            *error = "Stream ended before a matrix was found";
            return NULL;
        }

    if(c != '[')
        {
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

    if(c == EOF)
        {
            // ERROR: Stream ended where a numeric value
            //        was expected. 
            *error = "Stream ended where a numeric value "
                     "was expected";
            return NULL;
        }

    double *matrix = malloc(sizeof(matrix[0]) * 64);

    if(matrix == NULL)
        {
            // ERROR: Insufficient memory.
            *error = "Insufficient memory";
            return NULL;
        }
    
    int capacity = 64, size = 0, 
        w = -1, i = 0, j = 0;

    if(c != ']')
        while(1)
            {
                if(!isdigit(c))
                    {
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
                if(size == capacity)
                    {
                        int new_capacity = capacity * 2;

                        double *temp = realloc(matrix, sizeof(double) * new_capacity);

                        if(temp == NULL)
                            {
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

                if(errno)
                    {
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

                if(c == ']' || c == ',')
                    {
                        // The matrix's row just ended.

                        if(w == -1)
                            // This was the first row.
                            w = i;
                        else
                            {
                                // This wasn't the first row,
                                // so it's possible that it's
                                // length is different from the
                                // previous ones.
                                assert(w > -1);

                                if(i != w)
                                    {
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

                if(c == EOF)
                    {
                        // ERROR: Stream ended inside a matrix, where 
                        //        either ',', ']' or a numeric value was
                        //        expected.
                        *error = "Stream ended inside a matrix, where either "
                                 "',', ']' or a numeric value was expected";
                        return NULL;
                    }
            }

    if(size == 0)
        {
            free(matrix);
            *error = "Empty matrix";
            return NULL;
        }

    // If the internal fragmentation is too much,
    // return a dynamic memory region with the
    // exact size instead of the buffer used to
    // build the matrix.
    int fragm_threshold = 30; // (It's a percentage)
    
    if(100.0 * size/capacity < fragm_threshold)
        {
            int new_capacity = (size == 0) ? 1 : size;

            double *temp = realloc(matrix, new_capacity * sizeof(double));

            if(temp != NULL)
                matrix = temp;
        }

    *width = w;
    *height = j;

    return matrix;
}