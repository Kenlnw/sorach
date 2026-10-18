/*
 * sorach/matrix.h: 2-D matrices, part of the Sorach library.
 * https://github.com/Kenlnw/sorach
 * Released under the MIT License; see LICENSE.
 *
 * A SorMat is a rows x cols grid of SorScalar values stored in one flat block,
 * row by row (row-major order). Element (r, c) is at data[r * cols + c].
 *
 * Conventions shared by every function:
 *   - Results are written into a matrix the caller provides ("out").
 *     Functions do not allocate: only sor_mat_new takes memory from an arena.
 *     The functions that take an arena besides sor_mat_new borrow scratch
 *     space from it and give it back before they return.
 *   - A SorMat is passed by value. The struct holds only a shape and a
 *     pointer, so a copy still refers to the same elements.
 *   - Caller mistakes (wrong shape, NULL data, index out of range) stop the
 *     program with a message naming the failed check. These checks stay on
 *     in every build, including with -DNDEBUG. Only the file functions and
 *     the functions that can meet a singular matrix report failure through
 *     a return value: 1 for success, 0 for failure.
 *   - Row and column indexes start at 0.
 *
 * Typical use:
 *     SorArena arena = sor_arena_create(1024 * 1024);
 *
 *     SorMat a = sor_mat_new(&arena, 2, 3);
 *     SorMat b = sor_mat_new(&arena, 3, 2);
 *     SorMat c = sor_mat_new(&arena, 2, 2);
 *     sor_mat_fillall(a, 1);
 *     sor_mat_fillall(b, 2);
 *     sor_mat_mul(c, a, b);
 *     sor_mat_print(c);
 *
 *     sor_arena_destroy(&arena);
 */

#ifndef SORACH_MATRIX_H
#define SORACH_MATRIX_H

#include <stddef.h>
#include <stdio.h>
#include "core.h"
#include "arena.h"

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

typedef struct {
    int        rows, cols;   // shape
    SorScalar *data;         // rows * cols elements, row by row
} SorMat;

// Direction for the functions that work along one dimension. The axis names
// the dimension that is collapsed, as in NumPy, and the values are NumPy's
// axis numbers for a 2-D array:
//   SOR_MAT_ROW (axis 0): go down through the rows      -> one result per column
//   SOR_MAT_COL (axis 1): go across through the columns -> one result per row
typedef enum { SOR_MAT_ROW = 0, SOR_MAT_COL = 1 } SorMatAxis;

// Element (r, c) of matrix m, usable for both reading and writing:
//     SorScalar x = SOR_MAT_AT(m, 0, 2);
//     SOR_MAT_AT(m, 0, 2) = 5;
// No bounds checking is done. m is evaluated twice, so do not pass an
// expression with side effects.
#define SOR_MAT_AT(m, r, c) ((m).data[(r) * (m).cols + (c)])

// ---------------------------------------------------------------------------
// Creation
// ---------------------------------------------------------------------------

// Create a rows x cols matrix in the arena, with every element set to 0.
// rows and cols must be greater than 0. The matrix stays valid until the
// arena is reset, rewound past it, or destroyed.
SorMat sor_mat_new(SorArena *arena, int rows, int cols);

// ---------------------------------------------------------------------------
// Filling (set the contents of a matrix that already exists)
// ---------------------------------------------------------------------------

// Set every element to value.
void sor_mat_fillall(SorMat m, SorScalar value);

// Copy every element of src into dst. The shapes must match.
void sor_mat_copy(SorMat dst, SorMat src);

// Make out the identity matrix: 1 on the diagonal, 0 elsewhere.
// out must be square.
void sor_mat_identity(SorMat out);

// Copy rows * cols values from a plain array into out, row by row.
// src must hold at least that many values; this cannot be checked.
void sor_mat_from_array(SorMat out, const SorScalar *src);

// ---------------------------------------------------------------------------
// Random values
// ---------------------------------------------------------------------------

// Seed the random generator used by the functions below. The same seed
// always produces the same sequence, which makes runs repeatable.
// These functions use the C library's rand(), so they share its state with
// the rest of the program: calling sor_mat_seed is the same as calling srand,
// and any other code that calls rand() or srand() changes the sequence the
// functions below produce.
void sor_mat_seed(unsigned int seed);

// Fill out with values spread evenly over the range low <= x < high.
// low must be less than high.
void sor_mat_rand_uniform(SorMat out, SorScalar low, SorScalar high);

// Fill out with values from a normal (bell curve) distribution with the
// given mean and standard deviation. std must not be negative.
void sor_mat_rand_normal(SorMat out, SorScalar mean, SorScalar std);

// ---------------------------------------------------------------------------
// Views
// A view is a SorMat that points into another matrix's memory. Nothing is
// copied: writing through the view changes the original, and the view
// becomes invalid when the original's arena is reset, rewound, or destroyed.
// ---------------------------------------------------------------------------

// A 1 x cols view of one row of m.
SorMat sor_mat_row_view(SorMat m, int row);

// A count x cols view of consecutive rows of m, beginning at row start.
// The range start .. start + count - 1 must lie inside m.
SorMat sor_mat_rows_view(SorMat m, int start, int count);

// A view of all of m's elements with a different shape.
// rows * cols must equal m.rows * m.cols.
SorMat sor_mat_reshape(SorMat m, int rows, int cols);

// ---------------------------------------------------------------------------
// Rows and columns
// The get functions copy: changing out afterward does not change m.
// ---------------------------------------------------------------------------

// Copy column col of m into out. out must be m.rows x 1.
// (A column cannot be a view, because its elements are not adjacent.)
void sor_mat_get_col(SorMat out, SorMat m, int col);

// Copy src into column col of m. src must be m.rows x 1.
void sor_mat_set_col(SorMat m, int col, SorMat src);

// Copy row of m into out. out must be 1 x m.cols.
// To work on a row without copying it, use sor_mat_row_view.
void sor_mat_get_row(SorMat out, SorMat m, int row);

// Copy src into row of m. src must be 1 x m.cols.
void sor_mat_set_row(SorMat m, int row, SorMat src);

// Exchange two rows of m. Swapping a row with itself does nothing.
void sor_mat_swap_rows(SorMat m, int row_a, int row_b);

// Exchange two columns of m. Swapping a column with itself does nothing.
void sor_mat_swap_cols(SorMat m, int col_a, int col_b);

// Return 1 if a and b have the same number of rows and columns, else 0.
int sor_mat_same_shape(SorMat a, SorMat b);

// ---------------------------------------------------------------------------
// Linear algebra
// out must NOT share memory with an input: these functions read elements
// after writing others, so overlapping would corrupt the result.
// ---------------------------------------------------------------------------

// Matrix product: out = a @ b.
// Requires a.cols == b.rows; out must be a.rows x b.cols.
void sor_mat_mul(SorMat out, SorMat a, SorMat b);

// Matrix product with the first matrix transposed: out = a.T @ b.
// Gives the same result as sor_mat_transpose followed by sor_mat_mul, without
// needing a matrix to hold a.T.
// Requires a.rows == b.rows; out must be a.cols x b.cols.
void sor_mat_mul_tp(SorMat out, SorMat a, SorMat b);

// Transpose: out = a.T. out must be a.cols x a.rows.
void sor_mat_transpose(SorMat out, SorMat a);

// Transpose a square matrix in place.
void sor_mat_transpose_inplace(SorMat m);

// Make out a diagonal matrix with the elements of v on its diagonal and 0
// everywhere else. v may be a row (1 x n) or a column (n x 1) vector;
// out must be n x n.
void sor_mat_diag(SorMat out, SorMat v);

// Sum of the products of corresponding elements. a and b must have the same
// shape. For two vectors this is the dot product.
SorScalar sor_mat_dot(SorMat a, SorMat b);

// Sum of the diagonal elements. m must be square.
SorScalar sor_mat_trace(SorMat m);

// ---------------------------------------------------------------------------
// Elementwise
// All matrices involved must have the same shape. out MAY be the same
// matrix as an input, so sor_mat_add(a, a, b) computes a += b.
// ---------------------------------------------------------------------------

void sor_mat_add(SorMat out, SorMat a, SorMat b);        // out = a + b
void sor_mat_sub(SorMat out, SorMat a, SorMat b);        // out = a - b
void sor_mat_hadamard(SorMat out, SorMat a, SorMat b);   // out = a * b

// out = a / b. Division by a zero element is not checked and produces
// inf or nan, following normal floating-point rules.
void sor_mat_div(SorMat out, SorMat a, SorMat b);

// out = f(a): call f on every element of a. f must not be NULL.
void sor_mat_apply(SorMat out, SorMat a, SorScalar (*f)(SorScalar));

void sor_mat_scale(SorMat m, SorScalar s);                  // m *= s
void sor_mat_add_scalar(SorMat m, SorScalar s);             // m += s
void sor_mat_add_scaled(SorMat a, SorScalar s, SorMat b);   // a += s * b

// Limit every element to the range low .. high. low must not exceed high.
void sor_mat_clamp(SorMat m, SorScalar low, SorScalar high);

// ---------------------------------------------------------------------------
// Broadcasting a vector across a matrix
// The direction is decided by the shape of v:
//   v is 1 x m.cols: applied to every row
//   v is m.rows x 1: applied to every column
// Any other shape stops the program.
// ---------------------------------------------------------------------------

void sor_mat_add_vec(SorMat m, SorMat v);                // m += v
void sor_mat_sub_vec(SorMat m, SorMat v);                // m -= v
void sor_mat_mul_vec(SorMat m, SorMat v);                // m *= v
void sor_mat_div_vec(SorMat m, SorMat v);                // m /= v

// ---------------------------------------------------------------------------
// Reductions along an axis
//   SOR_MAT_ROW: collapse the rows    -> out must be 1 x m.cols (one value per column)
//   SOR_MAT_COL: collapse the columns -> out must be m.rows x 1 (one value per row)
// out must not be m.
// ---------------------------------------------------------------------------

void sor_mat_sum(SorMat out, SorMat m, SorMatAxis axis);
void sor_mat_mean(SorMat out, SorMat m, SorMatAxis axis);

// Standard deviation, dividing by the number of values n (the population
// standard deviation, the same as NumPy's np.std). It is 0 when a row or
// column holds a single value.
void sor_mat_std(SorMat out, SorMat m, SorMatAxis axis);

void sor_mat_max(SorMat out, SorMat m, SorMatAxis axis);
void sor_mat_min(SorMat out, SorMat m, SorMatAxis axis);

// Position of the largest (argmax) or smallest (argmin) value along an axis.
//   SOR_MAT_ROW: out receives m.cols entries, the row index of the result
//                in each column
//   SOR_MAT_COL: out receives m.rows entries, the column index of the result
//                in each row
// out_len must be exactly that number of entries. When several values tie,
// the first position is reported.
void sor_mat_argmax(int *out, size_t out_len, SorMat m, SorMatAxis axis);
void sor_mat_argmin(int *out, size_t out_len, SorMat m, SorMatAxis axis);

// ---------------------------------------------------------------------------
// Reductions over the whole matrix
// ---------------------------------------------------------------------------

SorScalar sor_mat_sum_all(SorMat m);
SorScalar sor_mat_mean_all(SorMat m);
SorScalar sor_mat_max_all(SorMat m);
SorScalar sor_mat_min_all(SorMat m);

// Frobenius norm: the square root of the sum of the squared elements.
// For a vector this is its length.
SorScalar sor_mat_norm(SorMat m);

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

// Return 1 if a and b have the same shape and every pair of corresponding
// elements differs by less than epsilon; otherwise return 0. Use a small
// positive epsilon (about 1e-5 for float, 1e-9 for double) instead of testing
// for exact equality, because floating-point results carry rounding error.
int sor_mat_equal(SorMat a, SorMat b, SorScalar epsilon);

// ---------------------------------------------------------------------------
// Printing
// ---------------------------------------------------------------------------

// Print the elements to standard output in NumPy style:
//     [[   1.0000   2.0000 ]
//      [   3.0000   4.0000 ]]
void sor_mat_print(SorMat m);

// Print the shape as "(rows, cols)".
void sor_mat_print_shape(SorMat m);

// ---------------------------------------------------------------------------
// Saving and loading
//
// Every function here uses the same binary format. Each stored matrix is a
// header (the bytes "SMAT", a format version, rows, cols and the element
// size) followed by its elements, row by row. Loading checks every part of
// the header, so a file of the wrong kind, shape or scalar type is reported
// instead of read as garbage. Numbers are stored in the byte order of the
// machine that wrote the file; a file from a machine with the opposite byte
// order is detected and refused.
//
// Every function here reports failure by printing a message to stderr and
// returning 0; success returns 1.
//
// sor_mat_save and sor_mat_load open and close the file themselves. The
// stream functions work on a FILE you opened, so several matrices can be
// stored one after another in one file. Open it in binary mode ("wb" or
// "rb"): on Windows, text mode changes some bytes and corrupts the data.
// ---------------------------------------------------------------------------

// Write m to the file at path, replacing any file already there.
int sor_mat_save(SorMat m, const char *path);

// Read a matrix saved by sor_mat_save or sor_mat_write into out. out must
// already have the stored shape, and the file must have been written with
// the same SorScalar type. On failure, out may be partly overwritten.
int sor_mat_load(SorMat out, const char *path);

// Write m to an open binary stream, at its current position. Does not open
// or close the stream.
int sor_mat_write(SorMat m, FILE *fp);

// Read the next stored matrix from an open binary stream into out, which
// must already have the stored shape. On failure, out may be partly
// overwritten.
int sor_mat_read_into(SorMat out, FILE *fp);

// Read the next stored matrix from an open binary stream, creating it in the
// arena with the stored shape. On success, sets *out. On failure, leaves
// *out unchanged and gives back any arena memory it took, so the arena is
// exactly as it was before the call.
int sor_mat_read(SorArena *arena, SorMat *out, FILE *fp);

// ---------------------------------------------------------------------------
// Square systems
// These use Gaussian elimination with partial pivoting. They need working
// space, which they take from the arena and give back before returning, so
// the arena is left exactly as it was.
// ---------------------------------------------------------------------------

// Determinant of m. m must be square and is not modified.
// Returns 0 for a singular matrix.
SorScalar sor_mat_det(SorArena *arena, SorMat m);

// Inverse: out = m^-1. m and out must be square with the same shape, and
// out must not be m. Returns 1 on success, or 0 if m is singular (has no
// inverse), in which case out is left unspecified.
int sor_mat_inverse(SorArena *arena, SorMat out, SorMat m);

// Solve A @ x = b for x. A must be n x n; b and x must both be n x k, so
// several right-hand sides can be solved at once. A and b are not modified.
// Returns 1 on success, or 0 if A is singular.
int sor_mat_solve(SorArena *arena, SorMat x, SorMat A, SorMat b);

#ifdef __cplusplus
}
#endif

#endif