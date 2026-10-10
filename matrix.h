#ifndef MATRIX_H
#define MATRIX_H

#include "lib/arena.h"

#ifdef MAT_USE_DOUBLE
typedef double Scalar;
#else
typedef float Scalar;
#endif

typedef struct {
    int rows, cols;
    Scalar *data;
} Matrix;

typedef enum { ROW = 0, COL = 1 } Axis;

#define MAT_AT(m, r, c) ((m).data[(r) * (m).cols + (c)])

// ---------------------------------------------------------------------------
// Allocation (the only function that uses an arena)
// ---------------------------------------------------------------------------
Matrix mat_new(Arena *arena, int rows, int cols);    // zero-filled

// ---------------------------------------------------------------------------
// Construction (fill a matrix that already exists)
// ---------------------------------------------------------------------------
void mat_fillall(Matrix m, Scalar value);           // every element = value
void mat_copy(Matrix dst, Matrix src);               // dst = src
void mat_identity(Matrix out);                       // 1s on the diagonal, 0s elsewhere (square)
void mat_from_array(Matrix out, const Scalar *src);   // copy rows*cols values from a plain array
void mat_rand_uniform(Matrix out, Scalar low, Scalar high);
void mat_rand_normal(Matrix out, Scalar mean, Scalar std);

// ---------------------------------------------------------------------------
// Views (no copying; the result shares memory with the original)
// ---------------------------------------------------------------------------
Matrix mat_row_view(Matrix m, int row);               // 1×cols matrix pointing at one row
Matrix mat_rows_view(Matrix m, int start, int count); // several consecutive rows
Matrix mat_reshape(Matrix m, int rows, int cols);     // same data, new shape

// ---------------------------------------------------------------------------
// Linear algebra
// out must NOT be the same matrix as a or b
// ---------------------------------------------------------------------------
void mat_mul(Matrix out, Matrix a, Matrix b);        // out = a @ b
void mat_transpose(Matrix out, Matrix a);            // out = a.T

// ---------------------------------------------------------------------------
// Elementwise
// out MAY be the same matrix as a or b
// ---------------------------------------------------------------------------
void mat_add(Matrix out, Matrix a, Matrix b);        // out = a + b
void mat_sub(Matrix out, Matrix a, Matrix b);        // out = a - b
void mat_hadamard(Matrix out, Matrix a, Matrix b);   // out = a * b per element
void mat_div(Matrix out, Matrix a, Matrix b);        // out = a / b per element
void mat_apply(Matrix out, Matrix a, Scalar (*f)(Scalar));  // out = f(a) per element

void mat_scale(Matrix m, Scalar s);             // mat *= scale
void mat_add_scalar(Matrix m, Scalar s);            // mat += s
void mat_add_scaled(Matrix y, Scalar s, Matrix x);// y += s * x

// ---------------------------------------------------------------------------
// Broadcasting a vector across a matrix (direction decided by v's shape)
//   v is 1×cols: applied to every row
//   v is rows×1: applied to every column
// ---------------------------------------------------------------------------
void mat_add_vec(Matrix m, Matrix v);                // m += v
void mat_sub_vec(Matrix m, Matrix v);                // m -= v
void mat_mul_vec(Matrix m, Matrix v);                // m *= v
void mat_div_vec(Matrix m, Matrix v);                // m /= v

// ---------------------------------------------------------------------------
// Reductions along an axis
//   axis ROW: collapse the rows    -> out is 1×cols (one value per column)
//   axis COL: collapse the columns -> out is rows×1 (one value per row)
// ---------------------------------------------------------------------------
void mat_sum(Matrix out, Matrix m, Axis axis);
void mat_mean(Matrix out, Matrix m, Axis axis);
void mat_max(Matrix out, Matrix m, Axis axis);
void mat_min(Matrix out, Matrix m, Axis axis);

// Index of the largest value along an axis
//   axis ROW: out has m.cols entries (row index of the max in each column)
//   axis COL: out has m.rows entries (column index of the max in each row)
void mat_argmax(int *out, size_t out_len, Matrix m, Axis axis);
void mat_argmin(int *out, size_t out_len, Matrix m, Axis axis);

// ---------------------------------------------------------------------------
// Reductions over the whole matrix
// ---------------------------------------------------------------------------
Scalar mat_sum_all(Matrix m);
Scalar mat_mean_all(Matrix m);
Scalar mat_max_all(Matrix m);
Scalar mat_min_all(Matrix m);
Scalar mat_norm(Matrix m);                            // sqrt of the sum of squares

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------
int mat_equal(Matrix a, Matrix b, Scalar eps);        // 1 if all elements differ by < eps

// ---------------------------------------------------------------------------
// Input and output
// ---------------------------------------------------------------------------
void mat_print(Matrix m);
void mat_save(Matrix m, const char *path);
void mat_load(Matrix out, const char *path);

// ---------------------------------------------------------------------------
// Square-matrix linear algebra (not needed for neural networks)
// ---------------------------------------------------------------------------
Scalar mat_trace(Matrix m);
Scalar mat_det(Matrix m);
void  mat_inverse(Matrix out, Matrix m);
void  mat_solve(Matrix x, Matrix A, Matrix b);       // solve A·x = b

#endif