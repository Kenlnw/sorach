#ifndef SORMAT_MATRIX_H
#define SORMAT_MATRIX_H

#include <stddef.h>
#include "arena.h"

/*
 * Matrix library.
 *
 * A Matrix is a rows x cols grid of MatScalar values stored in one flat block,
 * row by row (row-major order). Element (r, c) is at data[r * cols + c].
 *
 * Conventions shared by every function:
 *   - Results are written into a matrix the caller provides ("out").
 *     Functions do not allocate: only mat_new takes memory from an arena.
 *     The square-matrix functions borrow scratch space from an arena and
 *     give it back before they return.
 *   - A Matrix is passed by value. The struct holds only a shape and a
 *     pointer, so a copy still refers to the same elements.
 *   - Caller mistakes (wrong shape, NULL data, index out of range) stop the
 *     program with an assert. Only the file functions and the functions that
 *     can meet a singular matrix report failure through a return value.
 *   - Row and column indexes start at 0.
 *
 * Typical use:
 *     Arena arena = arena_create(1024 * 1024);
 *
 *     Matrix a = mat_new(&arena, 2, 3);
 *     Matrix b = mat_new(&arena, 3, 2);
 *     Matrix c = mat_new(&arena, 2, 2);
 *     mat_fillall(a, 1);
 *     mat_fillall(b, 2);
 *     mat_mul(c, a, b);
 *     mat_print(c);
 *
 *     arena_destroy(&arena);
 */

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

// The element type. float by default; compile with -DMAT_USE_DOUBLE to make
// every matrix use double. All files in a program must use the same setting.
#ifdef MAT_USE_DOUBLE
typedef double MatScalar;
#else
typedef float MatScalar;
#endif

typedef struct {
    int     rows, cols;   // shape
    MatScalar *data;         // rows * cols elements, row by row
} Matrix;

// Direction for the functions that work along one dimension. The axis names
// the dimension that is collapsed, as in NumPy:
//   MAT_ROW: go down through the rows    -> one result per column
//   MAT_COL: go across through the columns -> one result per row
typedef enum { MAT_ROW = 0, MAT_COL = 1 } MatAxis;

// Element (r, c) of matrix m, usable for both reading and writing:
//     MatScalar x = MAT_AT(m, 0, 2);
//     MAT_AT(m, 0, 2) = 5;
// No bounds checking is done. m is evaluated twice, so do not pass an
// expression with side effects.
#define MAT_AT(m, r, c) ((m).data[(r) * (m).cols + (c)])

// ---------------------------------------------------------------------------
// Creation
// ---------------------------------------------------------------------------

// Create a rows x cols matrix in the arena, with every element set to 0.
// rows and cols must be greater than 0. The matrix stays valid until the
// arena is reset, rewound past it, or destroyed.
Matrix mat_new(Arena *arena, int rows, int cols);

// ---------------------------------------------------------------------------
// Filling (set the contents of a matrix that already exists)
// ---------------------------------------------------------------------------

// Set every element to value.
void mat_fillall(Matrix m, MatScalar value);

// Copy every element of src into dst. The shapes must match.
void mat_copy(Matrix dst, Matrix src);

// Make out the identity matrix: 1 on the diagonal, 0 elsewhere.
// out must be square.
void mat_identity(Matrix out);

// Copy rows * cols values from a plain array into out, row by row.
// src must hold at least that many values; this cannot be checked.
void mat_from_array(Matrix out, const MatScalar *src);

// ---------------------------------------------------------------------------
// Random values
// ---------------------------------------------------------------------------

// Seed the random generator used by the functions below. The same seed
// always produces the same sequence, which makes runs repeatable.
void mat_seed(unsigned int seed);

// Fill out with values spread evenly over the range low <= x < high.
// low must be less than high.
void mat_rand_uniform(Matrix out, MatScalar low, MatScalar high);

// Fill out with values from a normal (bell curve) distribution with the
// given mean and standard deviation. std must not be negative.
void mat_rand_normal(Matrix out, MatScalar mean, MatScalar std);

// ---------------------------------------------------------------------------
// Views
// A view is a Matrix that points into another matrix's memory. Nothing is
// copied: writing through the view changes the original, and the view
// becomes invalid when the original's arena is reset, rewound, or destroyed.
// ---------------------------------------------------------------------------

// A 1 x cols view of one row of m.
// (The return matrix uses the SAME memory with m matrix. Please be careful.)
Matrix mat_row_view(Matrix m, int row);

// A count x cols view of consecutive rows of m, beginning at row start.
// The range start .. start + count - 1 must lie inside m.
// (The return matrix uses the SAME memory with m matrix. Please be careful.)
Matrix mat_rows_view(Matrix m, int start, int count);

// A view of all of m's elements with a different shape.
// rows * cols must equal m.rows * m.cols.
Matrix mat_reshape(Matrix m, int rows, int cols);

// ---------------------------------------------------------------------------
// Rows and columns
// ---------------------------------------------------------------------------

// Copy column col of m into out. out must be m.rows x 1.
// (A column cannot be a view, because its elements are not adjacent.)
void mat_get_col(Matrix out, Matrix m, int col);

// Copy src into column col of m. src must be m.rows x 1.
void mat_set_col(Matrix m, int col, Matrix src);

// Copy row of m into out. out must be 1 x m.cols.
// (The out matrix will not use the same memory with m matrix. Changing out doesn't matter m)
void mat_get_row(Matrix out, Matrix m, int row);

// Copy src into row of m. src must be 1 x m.cols.
void mat_set_row(Matrix m, int row, Matrix src);

// Exchange two rows of m. Swapping a row with itself does nothing.
void mat_swap_rows(Matrix m, int row_a, int row_b);

// Exchange two columns of m. Swapping a column with itself does nothing.
void mat_swap_cols(Matrix m, int col_a, int col_b);

// Return 1 if a and b have the same number of rows and columns, else 0.
int mat_same_shape(Matrix a, Matrix b);

// ---------------------------------------------------------------------------
// Linear algebra
// out must NOT share memory with an input: these functions read elements
// after writing others, so overlapping would corrupt the result.
// ---------------------------------------------------------------------------

// Matrix product: out = a @ b.
// Requires a.cols == b.rows; out must be a.rows x b.cols.
void mat_mul(Matrix out, Matrix a, Matrix b);

// Transpose: out = a.T. out must be a.cols x a.rows.
void mat_transpose(Matrix out, Matrix a);

// Transpose a square matrix in place.
void mat_transpose_inplace(Matrix m);

// Build a diagonal matrix from a row or column vector
void mat_diag(Matrix out, Matrix v);

// Sum of the products of corresponding elements. a and b must have the same
// shape. For two vectors this is the dot product.
MatScalar mat_dot(Matrix a, Matrix b);

// ---------------------------------------------------------------------------
// Elementwise
// All matrices involved must have the same shape. out MAY be the same
// matrix as an input, so mat_add(a, a, b) computes a += b.
// ---------------------------------------------------------------------------

void mat_add(Matrix out, Matrix a, Matrix b);        // out = a + b
void mat_sub(Matrix out, Matrix a, Matrix b);        // out = a - b
void mat_hadamard(Matrix out, Matrix a, Matrix b);   // out = a * b

// out = a / b. Division by a zero element is not checked and produces
// inf or nan, following normal floating-point rules.
void mat_div(Matrix out, Matrix a, Matrix b);

// out = f(a): call f on every element of a. f must not be NULL.
void mat_apply(Matrix out, Matrix a, MatScalar (*f)(MatScalar));

void mat_scale(Matrix m, MatScalar s);                  // m *= s
void mat_add_scalar(Matrix m, MatScalar s);             // m += s
void mat_add_scaled(Matrix a, MatScalar s, Matrix b);   // a += s * b

// Limit every element to the range low .. high. low must not exceed high.
void mat_clamp(Matrix m, MatScalar low, MatScalar high);

// ---------------------------------------------------------------------------
// Broadcasting a vector across a matrix
// The direction is decided by the shape of v:
//   v is 1 x m.cols: applied to every row
//   v is m.rows x 1: applied to every column
// Any other shape stops the program.
// ---------------------------------------------------------------------------

void mat_add_vec(Matrix m, Matrix v);                // m += v
void mat_sub_vec(Matrix m, Matrix v);                // m -= v
void mat_mul_vec(Matrix m, Matrix v);                // m *= v
void mat_div_vec(Matrix m, Matrix v);                // m /= v

// ---------------------------------------------------------------------------
// Reductions along an axis
//   MAT_ROW: collapse the rows    -> out must be 1 x m.cols (one value per column)
//   MAT_COL: collapse the columns -> out must be m.rows x 1 (one value per row)
// out must not be m.
// ---------------------------------------------------------------------------

void mat_sum(Matrix out, Matrix m, MatAxis axis);
void mat_mean(Matrix out, Matrix m, MatAxis axis);
void mat_max(Matrix out, Matrix m, MatAxis axis);
void mat_min(Matrix out, Matrix m, MatAxis axis);

// Position of the largest (argmax) or smallest (argmin) value along an axis.
//   ROW: out receives m.cols entries, the row index of the result in each column
//   COL: out receives m.rows entries, the column index of the result in each row
// out_len must be exactly that number of entries. When several values tie,
// the first position is reported.
void mat_argmax(int *out, size_t out_len, Matrix m, MatAxis axis);
void mat_argmin(int *out, size_t out_len, Matrix m, MatAxis axis);

// ---------------------------------------------------------------------------
// Reductions over the whole matrix
// ---------------------------------------------------------------------------

MatScalar mat_sum_all(Matrix m);
MatScalar mat_mean_all(Matrix m);
MatScalar mat_max_all(Matrix m);
MatScalar mat_min_all(Matrix m);

// Frobenius norm: the square root of the sum of the squared elements.
// For a vector this is its length.
MatScalar mat_norm(Matrix m);

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

// Return 1 if a and b have the same shape and every pair of corresponding
// elements differs by less than epsilon; otherwise return 0. Use a small
// positive epsilon (about 1e-5 for float, 1e-9 for double) instead of testing
// for exact equality, because floating-point results carry rounding error.
int mat_equal(Matrix a, Matrix b, MatScalar epsilon);

// ---------------------------------------------------------------------------
// Input and output
// ---------------------------------------------------------------------------

// Print the elements to standard output in NumPy style:
//     [[   1.0000   2.0000 ]
//      [   3.0000   4.0000 ]]
void mat_print(Matrix m);

// Print the shape as "(rows, cols)".
void mat_print_shape(Matrix m);

// Write m to a binary file: rows, cols and the element size as ints, then
// the elements. Returns 1 on success. On failure, prints a message to stderr
// and returns 0.
int mat_save(Matrix m, const char *path);

// Read a file written by mat_save into out. out must already have the shape
// stored in the file, and the file must have been written with the same
// MatScalar type. Returns 1 on success. On failure, prints a message to stderr
// and returns 0; out may then be partly overwritten.
int mat_load(Matrix out, const char *path);

// ---------------------------------------------------------------------------
// Square-matrix linear algebra (not needed for neural networks)
// These use Gaussian elimination with partial pivoting. They need working
// space, which they take from the arena and release before
// returning, so the arena is left exactly as it was.
// ---------------------------------------------------------------------------

// Sum of the diagonal elements. m must be square.
MatScalar mat_trace(Matrix m);

// Determinant of m. m must be square and is not modified.
// Returns 0 for a singular matrix.
MatScalar mat_det(Arena *arena, Matrix m);

// Inverse: out = m^-1. m and out must be square with the same shape, and
// out must not be m. Returns 1 on success, or 0 if m is singular (has no
// inverse), in which case out is left unspecified.
int mat_inverse(Arena *arena, Matrix out, Matrix m);

// Solve A @ x = b for x. A must be n x n; b and x must both be n x k, so
// several right-hand sides can be solved at once. A and b are not modified.
// Returns 1 on success, or 0 if A is singular.
int mat_solve(Arena *arena, Matrix x, Matrix A, Matrix b);

#endif