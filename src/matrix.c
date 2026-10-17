#include <sorach/matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include "check.h"

#define MAT_PI 3.14159265358979323846

// uniform in [0, 1)
static double rand01(void) {
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

// standard normal (mean 0, std 1), by the Box-Muller method
static double randn(void) {
    double u1 = 1.0 - rand01();   // in (0, 1], so log(u1) is safe
    double u2 = rand01();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * MAT_PI * u2);
}

static int gauss_jordan(SorMat work, SorMat rhs) {
    for (int col = 0; col < work.cols; col++) {
        int pivot_row = col;
        // Find the row where the largest absolute value is inside
        // That row will become a pivot row for one column.
        for (int row = col + 1; row < work.rows; row++) {
            if (fabs(SOR_MAT_AT(work, row, col)) > fabs(SOR_MAT_AT(work, pivot_row, col)))
                pivot_row = row;
        }

        if (SOR_MAT_AT(work, pivot_row, col) == 0) return 0;
        // We need a pivot row == col (all pivots should be diagnal), for perform a gaussian elimination
        if (pivot_row != col) {
            // do a row swap for both
            sor_mat_swap_rows(work, col, pivot_row);
            sor_mat_swap_rows(rhs, col, pivot_row);
        }

        // Scale the pivot row of BOTH by 1/pivot
        SorScalar pivot = SOR_MAT_AT(work, col, col);
        sor_mat_scale(sor_mat_row_view(work, col), 1 / pivot);
        sor_mat_scale(sor_mat_row_view(rhs, col), 1 / pivot);

        // Do a elimination for BOTH
        for (int row = 0; row < work.rows; row++) {
            if (row == col) continue;

            SorScalar s = -SOR_MAT_AT(work, row, col);
            sor_mat_add_scaled(sor_mat_row_view(work, row), s, sor_mat_row_view(work, col));
            sor_mat_add_scaled(sor_mat_row_view(rhs, row), s, sor_mat_row_view(rhs, col));
        }
    }
    return 1;
}

SorMat sor_mat_new(SorArena *arena, int rows, int cols) {
    SOR_CHECK(arena->buff != NULL);
    SOR_CHECK(rows > 0);
    SOR_CHECK(cols > 0);

    size_t data_size = (size_t) rows * cols * sizeof(SorScalar);
    SorMat m = {
        .rows = rows,
        .cols = cols,
        .data = sor_arena_alloc_zero(arena, data_size, SOR_ARENA_DEFAULT_ALIGN)
    };

    SOR_CHECK(m.data != NULL);
    return m;
}

void sor_mat_fillall(SorMat m, SorScalar value) {
    SOR_CHECK(m.data != NULL);
    for (int i = 0; i < (m.rows * m.cols); i++)
        m.data[i] = value;
}

void sor_mat_copy(SorMat dst, SorMat src) {
    SOR_CHECK(src.data != NULL && dst.data != NULL);
    SOR_CHECK(dst.rows == src.rows);
    SOR_CHECK(dst.cols == src.cols);

    // copy all elements in src.data into dst.data
    memcpy(dst.data, src.data, (size_t) src.rows * src.cols * sizeof(SorScalar));
}

void sor_mat_identity(SorMat out) {
    SOR_CHECK(out.data != NULL);
    SOR_CHECK(out.rows == out.cols);

    sor_mat_fillall(out, 0.0f);
    for (int i = 0; i < out.rows; i++)
        SOR_MAT_AT(out, i, i) = 1.0f;
}

void sor_mat_from_array(SorMat out, const SorScalar *src) {
    SOR_CHECK(out.data != NULL && src != NULL);
    memcpy(out.data, src, (size_t) out.rows * out.cols * sizeof(SorScalar));
}

void sor_mat_seed(unsigned int seed) {
    srand(seed);
}

void sor_mat_rand_uniform(SorMat out, SorScalar low, SorScalar high) {
    SOR_CHECK(out.data != NULL && low < high);

    for (int i = 0; i < out.rows * out.cols; i++)
        out.data[i] = (SorScalar)(low + (high - low) * rand01());
}

void sor_mat_rand_normal(SorMat out, SorScalar mean, SorScalar std) {
    SOR_CHECK(out.data != NULL && std >= 0);

    for (int i = 0; i < out.rows * out.cols; i++)
        out.data[i] = (SorScalar)(mean + std * randn());
}

SorMat sor_mat_row_view(SorMat m, int row) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(row >= 0 && row < m.rows);

    SorMat m_view = {
        .rows = 1,
        .cols = m.cols,
        .data = m.data + (row * m.cols)
    };
    return m_view;
}

SorMat sor_mat_rows_view(SorMat m, int start, int count) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(start >= 0 && start < m.rows);
    SOR_CHECK(count > 0 && start + count <= m.rows);

    SorMat m_view = {
        .rows = count,
        .cols = m.cols,
        .data = m.data + (start * m.cols)
    };
    return m_view;
}

SorMat sor_mat_reshape(SorMat m, int rows, int cols) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(rows > 0 && cols > 0);
    SOR_CHECK(rows * cols == m.rows * m.cols);

    SorMat new_shape = {
        .rows = rows,
        .cols = cols,
        .data = m.data
    };
    return new_shape;
}

void sor_mat_get_col(SorMat out, SorMat m, int col) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(col >= 0 && col < m.cols);
    SOR_CHECK(out.rows == m.rows && out.cols == 1);

    for (int row = 0; row < m.rows; row++)
        SOR_MAT_AT(out, row, 0) = SOR_MAT_AT(m, row, col);
}

void sor_mat_set_col(SorMat m, int col, SorMat src) {
    SOR_CHECK(m.data != NULL && src.data != NULL);
    SOR_CHECK(col >= 0 && col < m.cols);
    SOR_CHECK(src.rows == m.rows && src.cols == 1);

    for (int row = 0; row < m.rows; row++)
        SOR_MAT_AT(m, row, col) = SOR_MAT_AT(src, row, 0);
}

void sor_mat_get_row(SorMat out, SorMat m, int row) {
    sor_mat_copy(out, sor_mat_row_view(m, row));
}

void sor_mat_set_row(SorMat m, int row, SorMat src) {
    sor_mat_copy(sor_mat_row_view(m, row), src);
}

void sor_mat_swap_rows(SorMat m, int row_a, int row_b) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(row_a >= 0 && row_a < m.rows);
    SOR_CHECK(row_b >= 0 && row_b < m.rows);

    if (row_a == row_b) return;

    for (int col = 0; col < m.cols; col++) {
        SorScalar temp = SOR_MAT_AT(m, row_a, col);
        SOR_MAT_AT(m, row_a, col) = SOR_MAT_AT(m, row_b, col);
        SOR_MAT_AT(m, row_b, col) = temp;
    }
}

void sor_mat_swap_cols(SorMat m, int col_a, int col_b) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(col_a >= 0 && col_a < m.cols);
    SOR_CHECK(col_b >= 0 && col_b < m.cols);

    if (col_a == col_b) return;

    for (int row = 0; row < m.rows; row++) {
        SorScalar temp = SOR_MAT_AT(m, row, col_a);
        SOR_MAT_AT(m, row, col_a) = SOR_MAT_AT(m, row, col_b);
        SOR_MAT_AT(m, row, col_b) = temp;
    }
}

int sor_mat_same_shape(SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL);
    return a.rows == b.rows && a.cols == b.cols;
}

void sor_mat_mul(SorMat out, SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != a.data && out.data != b.data);
    SOR_CHECK(a.cols == b.rows);
    SOR_CHECK(out.rows == a.rows && out.cols == b.cols);

    sor_mat_fillall(out, 0.0);
    for (int a_row = 0; a_row < a.rows; a_row++) {
        for (int a_col = 0; a_col < a.cols; a_col++) {
            SorScalar a_ij = SOR_MAT_AT(a, a_row, a_col);
            for (int b_col = 0; b_col < b.cols; b_col++)
                SOR_MAT_AT(out, a_row, b_col) += a_ij * SOR_MAT_AT(b, a_col, b_col);
        }
    }
}

void sor_mat_mul_tp(SorMat out, SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != a.data && out.data != b.data);
    SOR_CHECK(a.rows == b.rows);
    SOR_CHECK(out.rows == a.cols && out.cols == b.cols);

    sor_mat_fillall(out, 0.0);
    for (int a_row = 0; a_row < a.rows; a_row++) {
        for (int a_col = 0; a_col < a.cols; a_col++) {
            SorScalar a_ij = SOR_MAT_AT(a, a_row, a_col);
            for (int b_col = 0; b_col < b.cols; b_col++)
                SOR_MAT_AT(out, a_col, b_col) += a_ij * SOR_MAT_AT(b, a_row, b_col);
        }
    }
}

void sor_mat_transpose(SorMat out, SorMat a) {
    SOR_CHECK(a.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != a.data);
    SOR_CHECK(out.rows == a.cols && out.cols == a.rows);

    for (int row = 0; row < out.rows; row++)
        for (int col = 0; col < out.cols; col++)
            SOR_MAT_AT(out, row, col) = SOR_MAT_AT(a, col, row);
}

void sor_mat_transpose_inplace(SorMat m) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(m.rows == m.cols);

    for (int row = 0; row < m.rows; row++) {
        for (int col = row + 1; col < m.rows; col++) {
            SorScalar temp = SOR_MAT_AT(m, row, col);
            SOR_MAT_AT(m, row, col) = SOR_MAT_AT(m, col, row);
            SOR_MAT_AT(m, col, row) = temp;
        }
    }
}

void sor_mat_diag(SorMat out, SorMat v) {
    SOR_CHECK(out.data != NULL && v.data != NULL);
    SOR_CHECK((v.cols == 1 && v.rows == out.rows) || (v.rows == 1 && v.cols == out.rows));
    SOR_CHECK(out.rows == out.cols);

    sor_mat_fillall(out, 0.0);
    for (int i = 0; i < v.rows; i++) SOR_MAT_AT(out, i, i) = v.data[i];
}

SorScalar sor_mat_dot(SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL);
    SOR_CHECK(sor_mat_same_shape(a, b));

    SorScalar dot = 0;
    for (int i = 0; i < a.rows * a.cols; i++)
        dot += a.data[i] * b.data[i];
    return dot;
}

void sor_mat_add(SorMat out, SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL && out.data != NULL);
    SOR_CHECK(a.rows == b.rows && a.cols == b.cols);
    SOR_CHECK(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] + b.data[i];
}

void sor_mat_sub(SorMat out, SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL && out.data != NULL);
    SOR_CHECK(a.rows == b.rows && a.cols == b.cols);
    SOR_CHECK(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] - b.data[i];
}

void sor_mat_hadamard(SorMat out, SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL && out.data != NULL);
    SOR_CHECK(a.rows == b.rows && a.cols == b.cols);
    SOR_CHECK(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] * b.data[i];
}

void sor_mat_div(SorMat out, SorMat a, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL && out.data != NULL);
    SOR_CHECK(a.rows == b.rows && a.cols == b.cols);
    SOR_CHECK(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] / b.data[i];
}

void sor_mat_apply(SorMat out, SorMat a, SorScalar (*func)(SorScalar)) {
    SOR_CHECK(a.data != NULL && out.data != NULL && func != NULL);
    SOR_CHECK(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = func(a.data[i]);
}

void sor_mat_scale(SorMat m, SorScalar s) {
    SOR_CHECK(m.data != NULL);

    for (int i = 0; i < m.rows * m.cols; i++)
        m.data[i] *= s;
}

void sor_mat_add_scalar(SorMat m, SorScalar s) {
    SOR_CHECK(m.data != NULL);

    for (int i = 0; i < m.rows * m.cols; i++)
        m.data[i] += s;
}

void sor_mat_add_scaled(SorMat a, SorScalar s, SorMat b) {
    SOR_CHECK(a.data != NULL && b.data != NULL);
    SOR_CHECK(a.rows == b.rows && a.cols == b.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        a.data[i] += s * b.data[i];
}

void sor_mat_clamp(SorMat m, SorScalar low, SorScalar high) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(low <= high);

    for (int i = 0; i < m.rows * m.cols; i++) {
        if (m.data[i] > high) m.data[i] = high;
        else if (m.data[i] < low) m.data[i] = low;
    }
}

void sor_mat_add_vec(SorMat m, SorMat v) {
    SOR_CHECK(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) += SOR_MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) += SOR_MAT_AT(v, m_row, 0);
        }
    }
    else SOR_CHECK(0 && "sor_mat_add_vec: v's size must be 1 x cols or rows x 1");
}

void sor_mat_sub_vec(SorMat m, SorMat v) {
    SOR_CHECK(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) -= SOR_MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) -= SOR_MAT_AT(v, m_row, 0);
        }
    }
    else SOR_CHECK(0 && "sor_mat_sub_vec: v's size must be 1 x cols or rows x 1");
}

void sor_mat_mul_vec(SorMat m, SorMat v) {
    SOR_CHECK(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) *= SOR_MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) *= SOR_MAT_AT(v, m_row, 0);
        }
    }
    else SOR_CHECK(0 && "sor_mat_mul_vec: v's size must be 1 x cols or rows x 1");
}

void sor_mat_div_vec(SorMat m, SorMat v) {
    SOR_CHECK(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) /= SOR_MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                SOR_MAT_AT(m, m_row, m_col) /= SOR_MAT_AT(v, m_row, 0);
        }
    }
    else SOR_CHECK(0 && "sor_mat_div_vec: v's size must be 1 x cols or rows x 1");
}

void sor_mat_sum(SorMat out, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != m.data);

    if (axis == SOR_MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            double sum = 0;
            for (int m_row = 0; m_row < m.rows; m_row++)
                sum += SOR_MAT_AT(m, m_row, m_col);
            SOR_MAT_AT(out, 0, m_col) = (SorScalar)sum;
        }
    }
    else if (axis == SOR_MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            double sum = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                sum += SOR_MAT_AT(m, m_row, m_col);
            SOR_MAT_AT(out, m_row, 0) = (SorScalar)sum;
        }
    }
    else SOR_CHECK(0 && "sor_mat_sum: out's size must be 1 x cols (SOR_MAT_ROW) or rows x 1 (SOR_MAT_COL)");
}

void sor_mat_mean(SorMat out, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != m.data);

    if (axis == SOR_MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            double sum = 0;
            for (int m_row = 0; m_row < m.rows; m_row++)
                sum += SOR_MAT_AT(m, m_row, m_col);
            SOR_MAT_AT(out, 0, m_col) = (SorScalar)(sum / m.rows);
        }
    }
    else if (axis == SOR_MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            double sum = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                sum += SOR_MAT_AT(m, m_row, m_col);
            SOR_MAT_AT(out, m_row, 0) = (SorScalar)(sum / m.cols);
        }
    }
    else SOR_CHECK(0 && "sor_mat_mean: out's size must be 1 x cols (SOR_MAT_ROW) or rows x 1 (SOR_MAT_COL)");
}

void sor_mat_std(SorMat out, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != m.data);

    if (axis == SOR_MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            double mean = 0;
             for (int m_row = 0; m_row < m.rows; m_row++)
                mean += SOR_MAT_AT(m, m_row, m_col);
            mean /= m.rows;

            double std = 0;
            for (int m_row = 0; m_row < m.rows; m_row++) {
                double diff = SOR_MAT_AT(m, m_row, m_col) - mean;
                std += diff * diff;
            }
            SOR_MAT_AT(out, 0, m_col) = (SorScalar)sqrt(std / m.rows);
        }
    }
    else if (axis == SOR_MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            double mean = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                mean += SOR_MAT_AT(m, m_row, m_col);
            mean /= m.cols;

            double std = 0;
            for (int m_col = 0; m_col < m.rows; m_col++) {
                double diff = SOR_MAT_AT(m, m_row, m_col) - mean;
                std += diff * diff;
            }
            SOR_MAT_AT(out, m_row, 0) = (SorScalar)sqrt(std / m.cols);
        }
    }
    else SOR_CHECK(0 && "sor_mat_std: out's size must be 1 x cols (SOR_MAT_ROW) or rows x 1 (SOR_MAT_COL)");
}

void sor_mat_max(SorMat out, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != m.data);

    if (axis == SOR_MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            SorScalar max = SOR_MAT_AT(m, 0, m_col);
            for (int m_row = 1; m_row < m.rows; m_row++)
                max = max < SOR_MAT_AT(m, m_row, m_col) ? SOR_MAT_AT(m, m_row, m_col) : max;
            SOR_MAT_AT(out, 0, m_col) = max;
        }
    }
    else if (axis == SOR_MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            SorScalar max = SOR_MAT_AT(m, m_row, 0);
            for (int m_col = 1; m_col < m.cols; m_col++)
                max = max < SOR_MAT_AT(m, m_row, m_col) ? SOR_MAT_AT(m, m_row, m_col) : max;
            SOR_MAT_AT(out, m_row, 0) = max;
        }
    }
    else SOR_CHECK(0 && "sor_mat_max: out's size must be 1 x cols (SOR_MAT_ROW) or rows x 1 (SOR_MAT_COL)");
}

void sor_mat_min(SorMat out, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(out.data != m.data);

    if (axis == SOR_MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            SorScalar min = SOR_MAT_AT(m, 0, m_col);
            for (int m_row = 1; m_row < m.rows; m_row++)
                min = min > SOR_MAT_AT(m, m_row, m_col) ? SOR_MAT_AT(m, m_row, m_col) : min;
            SOR_MAT_AT(out, 0, m_col) = min;
        }
    }
    else if (axis == SOR_MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            SorScalar min = SOR_MAT_AT(m, m_row, 0);
            for (int m_col = 1; m_col < m.cols; m_col++)
                min = min > SOR_MAT_AT(m, m_row, m_col) ? SOR_MAT_AT(m, m_row, m_col) : min;
            SOR_MAT_AT(out, m_row, 0) = min;
        }
    }
    else SOR_CHECK(0 && "sor_mat_min: out's size must be 1 x cols (SOR_MAT_ROW) or rows x 1 (SOR_MAT_COL)");
}

void sor_mat_argmax(int *out, size_t out_len, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out != NULL);

    if (axis == SOR_MAT_ROW && out_len == (size_t) m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            int max_idx = 0;
            for (int m_row = 1; m_row < m.rows; m_row++)
                max_idx = SOR_MAT_AT(m, max_idx, m_col) < SOR_MAT_AT(m, m_row, m_col) ? m_row : max_idx;
            out[m_col] = max_idx;
        }
    }
    else if (axis == SOR_MAT_COL && out_len == (size_t) m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            int max_idx = 0;
            for (int m_col = 1; m_col < m.cols; m_col++)
                max_idx = SOR_MAT_AT(m, m_row, max_idx) < SOR_MAT_AT(m, m_row, m_col) ? m_col : max_idx;
            out[m_row] = max_idx;
        }
    }
    else SOR_CHECK(0 && "sor_mat_argmax: out's length must be m.cols (SOR_MAT_ROW) or m.rows (SOR_MAT_COL)");
}

void sor_mat_argmin(int *out, size_t out_len, SorMat m, SorMatAxis axis) {
    SOR_CHECK(m.data != NULL && out != NULL);

    if (axis == SOR_MAT_ROW && out_len == (size_t) m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            int min_idx = 0;
            for (int m_row = 1; m_row < m.rows; m_row++)
                min_idx = SOR_MAT_AT(m, min_idx, m_col) > SOR_MAT_AT(m, m_row, m_col) ? m_row : min_idx;
            out[m_col] = min_idx;
        }
    }
    else if (axis == SOR_MAT_COL && out_len == (size_t) m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            int min_idx = 0;
            for (int m_col = 1; m_col < m.cols; m_col++)
                min_idx = SOR_MAT_AT(m, m_row, min_idx) > SOR_MAT_AT(m, m_row, m_col) ? m_col : min_idx;
            out[m_row] = min_idx;
        }
    }
    else SOR_CHECK(0 && "sor_mat_argmin: out's length must be m.cols (SOR_MAT_ROW) or m.rows (SOR_MAT_COL)");
}

SorScalar sor_mat_sum_all(SorMat m)  {
    SOR_CHECK(m.data != NULL);
    double sum = 0;
    for (int i = 0; i < m.rows * m.cols; i++) sum += m.data[i];
    return sum;
}

SorScalar sor_mat_mean_all(SorMat m)  {
    SOR_CHECK(m.data != NULL);
    double sum = 0;
    int m_len = m.rows * m.cols;
    for (int i = 0; i < m_len; i++) sum += m.data[i];
    return (SorScalar)(sum / m_len);
}

SorScalar sor_mat_max_all(SorMat m) {
    SOR_CHECK(m.data != NULL);
    SorScalar max = m.data[0];
    for (int i = 1; i < m.rows * m.cols; i++)
        max = max < m.data[i] ? m.data[i] : max;
    return max;
}

SorScalar sor_mat_min_all(SorMat m) {
    SOR_CHECK(m.data != NULL);
    SorScalar min = m.data[0];
    for (int i = 1; i < m.rows * m.cols; i++)
        min = min > m.data[i] ? m.data[i] : min;
    return min;
}

SorScalar sor_mat_norm(SorMat m) {
    SOR_CHECK(m.data != NULL);
    double squared_sum = 0;
    int m_len = m.rows * m.cols;
    for (int i = 0; i < m_len; i++)
        squared_sum += m.data[i] * m.data[i];
    return (SorScalar) sqrt(squared_sum);
}

int sor_mat_equal(SorMat a, SorMat b, SorScalar epsilon) {
    SOR_CHECK(a.data != NULL && b.data != NULL);
    if (a.rows != b.rows || a.cols != b.cols) return 0;

    for (int i = 0; i < a.rows * a.cols; i++)
        if (!(fabs(a.data[i] - b.data[i]) < epsilon)) return 0;
    return 1;
}

void sor_mat_print(SorMat m) {
    SOR_CHECK(m.data != NULL);

    for (int row = 0; row < m.rows; row++) {
        printf("%s", row == 0 ? "[[" : " [");
        for (int col = 0; col < m.cols; col++)
            printf(" %8.4f", (double)SOR_MAT_AT(m, row, col));
        printf("%s\n", row == m.rows - 1 ? " ]]" : " ]");
    }
}

void sor_mat_print_shape(SorMat m) {
    SOR_CHECK(m.data != NULL);
    printf("(%d, %d)\n", m.rows, m.cols);
}

SorScalar sor_mat_trace(SorMat m) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(m.rows == m.cols);

    double trace = 0.0;
    for (int i = 0; i < m.rows; i++)
        trace += SOR_MAT_AT(m, i, i);
    return (SorScalar)trace;
}

SorScalar sor_mat_det(SorArena *arena, SorMat m) {
    SOR_CHECK(m.data != NULL);
    SOR_CHECK(m.rows == m.cols);

    size_t mark = sor_arena_mark(arena);

    SorMat work = sor_mat_new(arena, m.rows, m.cols);
    sor_mat_copy(work, m);

    double det = 1.0;
    for (int col = 0; col < work.cols; col++) {
        int pivot_row = col;
        // Find the row where the largest absolute value is inside
        // That row will become a pivot row for one column.
        for (int row = col + 1; row < work.rows; row++) {
            if (fabs(SOR_MAT_AT(work, row, col)) > fabs(SOR_MAT_AT(work, pivot_row, col)))
                pivot_row = row;
        }

        if (SOR_MAT_AT(work, pivot_row, col) == 0) {
            det = 0.0;  // That's mean this matrix is singular
            break;
        }

        // We need a pivot row == col (all pivots should be diagnal), for perform a gaussian elimination
        if (pivot_row != col) {
            sor_mat_swap_rows(work, col, pivot_row);
            det = -det; // each row swap flips the sign of the determinant
        }

        for (int row = col + 1; row < work.rows; row++) {
            SorScalar s = -(SOR_MAT_AT(work, row, col) / SOR_MAT_AT(work, col, col));
            sor_mat_add_scaled(sor_mat_row_view(work, row), s, sor_mat_row_view(work, col));
        }
        det *= SOR_MAT_AT(work, col, col);
    }
    sor_arena_rewind(arena, mark);
    return (SorScalar)det;
}

int sor_mat_inverse(SorArena *arena, SorMat out, SorMat m) {
    SOR_CHECK(m.data != NULL && out.data != NULL);
    SOR_CHECK(m.rows == m.cols);
    SOR_CHECK(sor_mat_same_shape(m, out));
    SOR_CHECK(out.data != m.data);

    size_t mark = sor_arena_mark(arena);

    SorMat work = sor_mat_new(arena, m.rows, m.cols);
    sor_mat_copy(work, m);
    sor_mat_identity(out);

    int is_completed = gauss_jordan(work, out);

    sor_arena_rewind(arena, mark);
    return is_completed;
}

int sor_mat_solve(SorArena *arena, SorMat x, SorMat A, SorMat b) {
    SOR_CHECK(x.data != NULL && A.data != NULL && b.data != NULL);
    SOR_CHECK(A.rows == A.cols);
    SOR_CHECK(A.rows == b.rows);
    SOR_CHECK(sor_mat_same_shape(b, x));
    SOR_CHECK(x.data != A.data && x.data != b.data);

    size_t mark = sor_arena_mark(arena);

    SorMat work = sor_mat_new(arena, A.rows, A.cols);
    sor_mat_copy(work, A);

    sor_mat_copy(x, b);

    int is_completed = gauss_jordan(work, x);

    sor_arena_rewind(arena, mark);
    return is_completed;
}