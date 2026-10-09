#include "matrix.h"
#include<stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

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

Matrix mat_new(Arena *arena, int rows, int cols) {
    assert(arena->buff != NULL);
    assert(rows > 0);
    assert(cols > 0);

    size_t data_size = (size_t) rows * cols * sizeof(Scalar);
    Matrix m = {
        .rows = rows,
        .cols = cols,
        .data = arena_alloc_zero(arena, data_size, DEFAULT_ALIGN)
    };

    assert(m.data != NULL);
    return m;
}

void mat_fillall(Matrix m, Scalar value) {
    assert(m.data != NULL);
    for (int i = 0; i < (m.rows * m.cols); i++)
        m.data[i] = value;
}

void mat_copy(Matrix dst, Matrix src) {
    assert(src.data != NULL && dst.data != NULL);
    assert(dst.rows == src.rows);
    assert(dst.cols == src.cols);

    // copy all elements in src.data into dst.data
    memcpy(dst.data, src.data, (size_t) src.rows * src.cols * sizeof(Scalar));
}

void mat_identity(Matrix out) {
    assert(out.data != NULL);
    assert(out.rows == out.cols);

    mat_fillall(out, 0.0f);
    for (int i = 0; i < out.rows; i++)
        MAT_AT(out, i, i) = 1.0f;
}

void mat_from_array(Matrix out, const Scalar *src) {
    assert(out.data != NULL && src != NULL);
    memcpy(out.data, src, (size_t) out.rows * out.cols * sizeof(Scalar));
}

void mat_rand_uniform(Matrix out, Scalar low, Scalar high) {
    assert(out.data != NULL && low < high);

    for (int i = 0; i < out.rows * out.cols; i++)
        out.data[i] = (Scalar)(low + (high - low) * rand01());
}

void mat_rand_normal(Matrix out, Scalar mean, Scalar std) {
    assert(out.data != NULL && std >= 0);

    for (int i = 0; i < out.rows * out.cols; i++)
        out.data[i] = (Scalar)(mean + std * randn());
}

Matrix mat_row_view(Matrix m, int row) {
    assert(m.data != NULL);
    assert(row >= 0 && row < m.rows);

    Matrix m_view = {
        .rows = 1,
        .cols = m.cols,
        .data = m.data + (row * m.cols)
    };
    return m_view;
}

Matrix mat_rows_view(Matrix m, int start, int count) {
    assert(m.data != NULL);
    assert(start >= 0 && start < m.rows);
    assert(count > 0 && start + count <= m.rows);

    Matrix m_view = {
        .rows = count,
        .cols = m.cols,
        .data = m.data + (start * m.cols)
    };
    return m_view;
}

Matrix mat_reshape(Matrix m, int rows, int cols) {
    assert(m.data != NULL);
    assert(rows > 0 && cols > 0);
    assert(rows * cols == m.rows * m.cols);

    Matrix new_shape = {
        .rows = rows,
        .cols = cols,
        .data = m.data
    };
    return new_shape;
}

void mat_mul(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(out.data != a.data && out.data != b.data);
    assert(a.cols == b.rows);
    assert(out.rows == a.rows && out.cols == b.cols);

    for (int a_row = 0; a_row < a.rows; a_row++) {
        for (int b_col = 0; b_col < b.cols; b_col++) {
            Scalar sum = 0;
            for (int a_col = 0; a_col < a.cols; a_col++)
                sum += MAT_AT(a, a_row, a_col) * MAT_AT(b, a_col, b_col);
            MAT_AT(out, a_row, b_col) = sum;
        }
    }
}

void mat_transpose(Matrix out, Matrix a) {
    assert(a.data != NULL && out.data != NULL);
    assert(out.data != a.data);
    assert(out.rows == a.cols && out.cols == a.rows);

    for (int row = 0; row < out.rows; row++)
        for (int col = 0; col < out.cols; col++)
            MAT_AT(out, row, col) = MAT_AT(a, col, row);
}

void mat_add(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(out.data != a.data && out.data != b.data);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] + b.data[i];
}

void mat_sub(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(out.data != a.data && out.data != b.data);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] - b.data[i];
}

void mat_hadamard(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(out.data != a.data && out.data != b.data);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] * b.data[i];
}

void mat_div(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(out.data != a.data && out.data != b.data);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] / b.data[i];
}

void mat_apply(Matrix out, Matrix a, Scalar (*func)(Scalar)) {
    assert(a.data != NULL && out.data != NULL && func != NULL);
    assert(out.data != a.data);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = func(a.data[i]);
}

void mat_scale(Matrix m, Scalar s) {
    assert(m.data != NULL);

    for (int i = 0; i < m.rows * m.cols; i++)
        m.data[i] *= s;
}

void mat_add_scalar(Matrix m, Scalar s) {
    assert(m.data != NULL);

    for (int i = 0; i < m.rows * m.cols; i++)
        m.data[i] += s;
}

void mat_add_scaled(Matrix a, Scalar s, Matrix b) {
    assert(a.data != NULL && b.data != NULL);
    assert(a.rows == b.rows && a.cols == b.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        a.data[i] += s * b.data[i];
}

void mat_add_vec(Matrix m, Matrix v) {
    assert(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                MAT_AT(m, m_row, m_col) += MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_col = 0; m_col < m.rows; m_col++) {
            for (int m_row = 0; m_row < m.cols; m_row++)
                MAT_AT(m, m_row, m_col) += MAT_AT(v, m_row, 0);
        }
    }
    else assert(0 && "mat_add_vec: v's size must be 1 x cols or rows x 1");
}

void mat_sub_vec(Matrix m, Matrix v) {
    assert(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                MAT_AT(m, m_row, m_col) -= MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_col = 0; m_col < m.rows; m_col++) {
            for (int m_row = 0; m_row < m.cols; m_row++)
                MAT_AT(m, m_row, m_col) -= MAT_AT(v, m_row, 0);
        }
    }
    else assert(0 && "mat_sub_vec: v's size must be 1 x cols or rows x 1");
}

void mat_mul_vec(Matrix m, Matrix v) {
    assert(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                MAT_AT(m, m_row, m_col) *= MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_col = 0; m_col < m.rows; m_col++) {
            for (int m_row = 0; m_row < m.cols; m_row++)
                MAT_AT(m, m_row, m_col) *= MAT_AT(v, m_row, 0);
        }
    }
    else assert(0 && "mat_mul_vec: v's size must be 1 x cols or rows x 1");
}

void mat_div_vec(Matrix m, Matrix v) {
    assert(m.data != NULL && v.data != NULL);

    if (v.rows == 1 && v.cols == m.cols) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                MAT_AT(m, m_row, m_col) /= MAT_AT(v, 0, m_col);
        }
    }
    else if (v.cols == 1 && v.rows == m.rows) {
        for (int m_col = 0; m_col < m.rows; m_col++) {
            for (int m_row = 0; m_row < m.cols; m_row++)
                MAT_AT(m, m_row, m_col) /= MAT_AT(v, m_row, 0);
        }
    }
    else assert(0 && "mat_div_vec: v's size must be 1 x cols or rows x 1");
}

void mat_sum(Matrix out, Matrix m, Axis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            Scalar sum = 0;
            for (int m_row = 0; m_row < m.rows; m_row++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, 0, m_col) = sum;
        }
    }
    else if (axis == COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            Scalar sum = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, m_row, 0) = sum;
        }
    }
    else assert(0 && "mat_sum: out's size must be 1 x cols (ROW) or rows x 1 (COL)");
}

void mat_mean(Matrix out, Matrix m, Axis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            Scalar sum = 0;
            for (int m_row = 0; m_row < m.rows; m_row++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, 0, m_col) = (Scalar)(sum / m.rows);
        }
    }
    else if (axis == COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            Scalar sum = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, m_row, 0) = (Scalar)(sum / m.cols);
        }
    }
    else assert(0 && "mat_mean: out's size must be 1 x cols (ROW) or rows x 1 (COL)");
}

void mat_max(Matrix out, Matrix m, Axis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            Scalar max = MAT_AT(m, 0, m_col);
            for (int m_row = 1; m_row < m.rows; m_row++)
                max = max < MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : max;
            MAT_AT(out, 0, m_col) = max;
        }
    }
    else if (axis == COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            Scalar max = MAT_AT(m, m_row, 0);
            for (int m_col = 1; m_col < m.cols; m_col++)
                max = max < MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : max;
            MAT_AT(out, m_row, 0) = max;
        }
    }
    else assert(0 && "mat_max: out's size must be 1 x cols (ROW) or rows x 1 (COL)");
}

void mat_min(Matrix out, Matrix m, Axis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            Scalar min = MAT_AT(m, 0, m_col);
            for (int m_row = 1; m_row < m.rows; m_row++)
                min = min > MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : min;
            MAT_AT(out, 0, m_col) = min;
        }
    }
    else if (axis == COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            Scalar min = MAT_AT(m, m_row, 0);
            for (int m_col = 1; m_col < m.cols; m_col++)
                min = min > MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : min;
            MAT_AT(out, m_row, 0) = min;
        }
    }
    else assert(0 && "mat_min: out's size must be 1 x cols (ROW) or rows x 1 (COL)");
}