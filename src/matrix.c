#include <sora/matrix.h>
#include <stdlib.h>
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

static int gauss_jordan(Matrix work, Matrix rhs) {
    for (int col = 0; col < work.cols; col++) {
        int pivot_row = col;
        // Find the row where the largest absolute value is inside
        // That row will become a pivot row for one column.
        for (int row = col + 1; row < work.rows; row++) {
            if (fabs(MAT_AT(work, row, col)) > fabs(MAT_AT(work, pivot_row, col)))
                pivot_row = row;
        }

        if (MAT_AT(work, pivot_row, col) == 0) return 0;
        // We need a pivot row == col (all pivots should be diagnal), for perform a gaussian elimination
        if (pivot_row != col) {
            // do a row swap for both
            mat_swap_rows(work, col, pivot_row);
            mat_swap_rows(rhs, col, pivot_row);
        }

        // Scale the pivot row of BOTH by 1/pivot
        MatScalar pivot = MAT_AT(work, col, col);
        mat_scale(mat_row_view(work, col), 1 / pivot);
        mat_scale(mat_row_view(rhs, col), 1 / pivot);

        // Do a elimination for BOTH
        for (int row = 0; row < work.rows; row++) {
            if (row == col) continue;

            MatScalar s = -MAT_AT(work, row, col);
            mat_add_scaled(mat_row_view(work, row), s, mat_row_view(work, col));
            mat_add_scaled(mat_row_view(rhs, row), s, mat_row_view(rhs, col));
        }
    }
    return 1;
}

Matrix mat_new(Arena *arena, int rows, int cols) {
    assert(arena->buff != NULL);
    assert(rows > 0);
    assert(cols > 0);

    size_t data_size = (size_t) rows * cols * sizeof(MatScalar);
    Matrix m = {
        .rows = rows,
        .cols = cols,
        .data = arena_alloc_zero(arena, data_size, ARENA_DEFAULT_ALIGN)
    };

    assert(m.data != NULL);
    return m;
}

void mat_fillall(Matrix m, MatScalar value) {
    assert(m.data != NULL);
    for (int i = 0; i < (m.rows * m.cols); i++)
        m.data[i] = value;
}

void mat_copy(Matrix dst, Matrix src) {
    assert(src.data != NULL && dst.data != NULL);
    assert(dst.rows == src.rows);
    assert(dst.cols == src.cols);

    // copy all elements in src.data into dst.data
    memcpy(dst.data, src.data, (size_t) src.rows * src.cols * sizeof(MatScalar));
}

void mat_identity(Matrix out) {
    assert(out.data != NULL);
    assert(out.rows == out.cols);

    mat_fillall(out, 0.0f);
    for (int i = 0; i < out.rows; i++)
        MAT_AT(out, i, i) = 1.0f;
}

void mat_from_array(Matrix out, const MatScalar *src) {
    assert(out.data != NULL && src != NULL);
    memcpy(out.data, src, (size_t) out.rows * out.cols * sizeof(MatScalar));
}

void mat_seed(unsigned int seed) {
    srand(seed);
}

void mat_rand_uniform(Matrix out, MatScalar low, MatScalar high) {
    assert(out.data != NULL && low < high);

    for (int i = 0; i < out.rows * out.cols; i++)
        out.data[i] = (MatScalar)(low + (high - low) * rand01());
}

void mat_rand_normal(Matrix out, MatScalar mean, MatScalar std) {
    assert(out.data != NULL && std >= 0);

    for (int i = 0; i < out.rows * out.cols; i++)
        out.data[i] = (MatScalar)(mean + std * randn());
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

void mat_get_col(Matrix out, Matrix m, int col) {
    assert(m.data != NULL && out.data != NULL);
    assert(col >= 0 && col < m.cols);
    assert(out.rows == m.rows && out.cols == 1);

    for (int row = 0; row < m.rows; row++)
        MAT_AT(out, row, 0) = MAT_AT(m, row, col);
}

void mat_set_col(Matrix m, int col, Matrix src) {
    assert(m.data != NULL && src.data != NULL);
    assert(col >= 0 && col < m.cols);
    assert(src.rows == m.rows && src.cols == 1);

    for (int row = 0; row < m.rows; row++)
        MAT_AT(m, row, col) = MAT_AT(src, row, 0);
}

void mat_get_row(Matrix out, Matrix m, int row) {
    mat_copy(out, mat_row_view(m, row));
}

void mat_set_row(Matrix m, int row, Matrix src) {
    mat_copy(mat_row_view(m, row), src);
}

void mat_swap_rows(Matrix m, int row_a, int row_b) {
    assert(m.data != NULL);
    assert(row_a >= 0 && row_a < m.rows);
    assert(row_b >= 0 && row_b < m.rows);

    if (row_a == row_b) return;

    for (int col = 0; col < m.cols; col++) {
        MatScalar temp = MAT_AT(m, row_a, col);
        MAT_AT(m, row_a, col) = MAT_AT(m, row_b, col);
        MAT_AT(m, row_b, col) = temp;
    }
}

void mat_swap_cols(Matrix m, int col_a, int col_b) {
    assert(m.data != NULL);
    assert(col_a >= 0 && col_a < m.cols);
    assert(col_b >= 0 && col_b < m.cols);

    if (col_a == col_b) return;

    for (int row = 0; row < m.rows; row++) {
        MatScalar temp = MAT_AT(m, row, col_a);
        MAT_AT(m, row, col_a) = MAT_AT(m, row, col_b);
        MAT_AT(m, row, col_b) = temp;
    }
}

int mat_same_shape(Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL);
    return a.rows == b.rows && a.cols == b.cols;
}

void mat_mul(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(out.data != a.data && out.data != b.data);
    assert(a.cols == b.rows);
    assert(out.rows == a.rows && out.cols == b.cols);

    for (int a_row = 0; a_row < a.rows; a_row++) {
        for (int b_col = 0; b_col < b.cols; b_col++) {
            MatScalar sum = 0;
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

void mat_transpose_inplace(Matrix m) {
    assert(m.data != NULL);
    assert(m.rows == m.cols);

    for (int row = 0; row < m.rows; row++) {
        for (int col = row + 1; col < m.rows; col++) {
            MatScalar temp = MAT_AT(m, row, col);
            MAT_AT(m, row, col) = MAT_AT(m, col, row);
            MAT_AT(m, col, row) = temp;
        }
    }
}

MatScalar mat_dot(Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL);
    assert(mat_same_shape(a, b));

    MatScalar dot = 0;
    for (int i = 0; i < a.rows * a.cols; i++)
        dot += a.data[i] * b.data[i];
    return dot;
}

void mat_add(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] + b.data[i];
}

void mat_sub(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] - b.data[i];
}

void mat_hadamard(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] * b.data[i];
}

void mat_div(Matrix out, Matrix a, Matrix b) {
    assert(a.data != NULL && b.data != NULL && out.data != NULL);
    assert(a.rows == b.rows && a.cols == b.cols);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = a.data[i] / b.data[i];
}

void mat_apply(Matrix out, Matrix a, MatScalar (*func)(MatScalar)) {
    assert(a.data != NULL && out.data != NULL && func != NULL);
    assert(out.rows == a.rows && out.cols == a.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        out.data[i] = func(a.data[i]);
}

void mat_scale(Matrix m, MatScalar s) {
    assert(m.data != NULL);

    for (int i = 0; i < m.rows * m.cols; i++)
        m.data[i] *= s;
}

void mat_add_scalar(Matrix m, MatScalar s) {
    assert(m.data != NULL);

    for (int i = 0; i < m.rows * m.cols; i++)
        m.data[i] += s;
}

void mat_add_scaled(Matrix a, MatScalar s, Matrix b) {
    assert(a.data != NULL && b.data != NULL);
    assert(a.rows == b.rows && a.cols == b.cols);

    for (int i = 0; i < a.rows * a.cols; i++)
        a.data[i] += s * b.data[i];
}

void mat_clamp(Matrix m, MatScalar low, MatScalar high) {
    assert(m.data != NULL);
    assert(low <= high);

    for (int i = 0; i < m.rows * m.cols; i++) {
        if (m.data[i] > high) m.data[i] = high;
        else if (m.data[i] < low) m.data[i] = low;
    }
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
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
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
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
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
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
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
        for (int m_row = 0; m_row < m.rows; m_row++) {
            for (int m_col = 0; m_col < m.cols; m_col++)
                MAT_AT(m, m_row, m_col) /= MAT_AT(v, m_row, 0);
        }
    }
    else assert(0 && "mat_div_vec: v's size must be 1 x cols or rows x 1");
}

void mat_sum(Matrix out, Matrix m, MatAxis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            double sum = 0;
            for (int m_row = 0; m_row < m.rows; m_row++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, 0, m_col) = (MatScalar)sum;
        }
    }
    else if (axis == MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            double sum = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, m_row, 0) = (MatScalar)sum;
        }
    }
    else assert(0 && "mat_sum: out's size must be 1 x cols (MAT_ROW) or rows x 1 (MAT_COL)");
}

void mat_mean(Matrix out, Matrix m, MatAxis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            double sum = 0;
            for (int m_row = 0; m_row < m.rows; m_row++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, 0, m_col) = (MatScalar)(sum / m.rows);
        }
    }
    else if (axis == MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            double sum = 0;
            for (int m_col = 0; m_col < m.cols; m_col++)
                sum += MAT_AT(m, m_row, m_col);
            MAT_AT(out, m_row, 0) = (MatScalar)(sum / m.cols);
        }
    }
    else assert(0 && "mat_mean: out's size must be 1 x cols (MAT_ROW) or rows x 1 (MAT_COL)");
}

void mat_max(Matrix out, Matrix m, MatAxis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            MatScalar max = MAT_AT(m, 0, m_col);
            for (int m_row = 1; m_row < m.rows; m_row++)
                max = max < MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : max;
            MAT_AT(out, 0, m_col) = max;
        }
    }
    else if (axis == MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            MatScalar max = MAT_AT(m, m_row, 0);
            for (int m_col = 1; m_col < m.cols; m_col++)
                max = max < MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : max;
            MAT_AT(out, m_row, 0) = max;
        }
    }
    else assert(0 && "mat_max: out's size must be 1 x cols (MAT_ROW) or rows x 1 (MAT_COL)");
}

void mat_min(Matrix out, Matrix m, MatAxis axis) {
    assert(m.data != NULL && out.data != NULL);
    assert(out.data != m.data);

    if (axis == MAT_ROW && out.rows == 1 && out.cols == m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            MatScalar min = MAT_AT(m, 0, m_col);
            for (int m_row = 1; m_row < m.rows; m_row++)
                min = min > MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : min;
            MAT_AT(out, 0, m_col) = min;
        }
    }
    else if (axis == MAT_COL && out.rows == m.rows && out.cols == 1) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            MatScalar min = MAT_AT(m, m_row, 0);
            for (int m_col = 1; m_col < m.cols; m_col++)
                min = min > MAT_AT(m, m_row, m_col) ? MAT_AT(m, m_row, m_col) : min;
            MAT_AT(out, m_row, 0) = min;
        }
    }
    else assert(0 && "mat_min: out's size must be 1 x cols (MAT_ROW) or rows x 1 (MAT_COL)");
}

void mat_argmax(int *out, size_t out_len, Matrix m, MatAxis axis) {
    assert(m.data != NULL && out != NULL);

    if (axis == MAT_ROW && out_len == (size_t) m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            int max_idx = 0;
            for (int m_row = 1; m_row < m.rows; m_row++)
                max_idx = MAT_AT(m, max_idx, m_col) < MAT_AT(m, m_row, m_col) ? m_row : max_idx;
            out[m_col] = max_idx;
        }
    }
    else if (axis == MAT_COL && out_len == (size_t) m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            int max_idx = 0;
            for (int m_col = 1; m_col < m.cols; m_col++)
                max_idx = MAT_AT(m, m_row, max_idx) < MAT_AT(m, m_row, m_col) ? m_col : max_idx;
            out[m_row] = max_idx;
        }
    }
    else assert(0 && "mat_argmax: out's length must be m.cols (MAT_ROW) or m.rows (MAT_COL)");
}

void mat_argmin(int *out, size_t out_len, Matrix m, MatAxis axis) {
    assert(m.data != NULL && out != NULL);

    if (axis == MAT_ROW && out_len == (size_t) m.cols) {
        for (int m_col = 0; m_col < m.cols; m_col++) {
            int min_idx = 0;
            for (int m_row = 1; m_row < m.rows; m_row++)
                min_idx = MAT_AT(m, min_idx, m_col) > MAT_AT(m, m_row, m_col) ? m_row : min_idx;
            out[m_col] = min_idx;
        }
    }
    else if (axis == MAT_COL && out_len == (size_t) m.rows) {
        for (int m_row = 0; m_row < m.rows; m_row++) {
            int min_idx = 0;
            for (int m_col = 1; m_col < m.cols; m_col++)
                min_idx = MAT_AT(m, m_row, min_idx) > MAT_AT(m, m_row, m_col) ? m_col : min_idx;
            out[m_row] = min_idx;
        }
    }
    else assert(0 && "mat_argmin: out's length must be m.cols (MAT_ROW) or m.rows (MAT_COL)");
}

MatScalar mat_sum_all(Matrix m)  {
    assert(m.data != NULL);
    double sum = 0;
    for (int i = 0; i < m.rows * m.cols; i++) sum += m.data[i];
    return sum;
}

MatScalar mat_mean_all(Matrix m)  {
    assert(m.data != NULL);
    double sum = 0;
    int m_len = m.rows * m.cols;
    for (int i = 0; i < m_len; i++) sum += m.data[i];
    return (MatScalar)(sum / m_len);
}

MatScalar mat_max_all(Matrix m) {
    assert(m.data != NULL);
    MatScalar max = m.data[0];
    for (int i = 1; i < m.rows * m.cols; i++)
        max = max < m.data[i] ? m.data[i] : max;
    return max;
}

MatScalar mat_min_all(Matrix m) {
    assert(m.data != NULL);
    MatScalar min = m.data[0];
    for (int i = 1; i < m.rows * m.cols; i++)
        min = min > m.data[i] ? m.data[i] : min;
    return min;
}

MatScalar mat_norm(Matrix m) {
    assert(m.data != NULL);
    double squared_sum = 0;
    int m_len = m.rows * m.cols;
    for (int i = 0; i < m_len; i++)
        squared_sum += m.data[i] * m.data[i];
    return (MatScalar) sqrt(squared_sum);
}

int mat_equal(Matrix a, Matrix b, MatScalar epsilon) {
    assert(a.data != NULL && b.data != NULL);
    if (a.rows != b.rows || a.cols != b.cols) return 0;

    for (int i = 0; i < a.rows * a.cols; i++)
        if (!(fabs(a.data[i] - b.data[i]) < epsilon)) return 0;
    return 1;
}

void mat_print(Matrix m) {
    assert(m.data != NULL);

    for (int row = 0; row < m.rows; row++) {
        printf("%s", row == 0 ? "[[" : " [");
        for (int col = 0; col < m.cols; col++)
            printf(" %8.4f", (double)MAT_AT(m, row, col));
        printf("%s\n", row == m.rows - 1 ? " ]]" : " ]");
    }
}

void mat_print_shape(Matrix m) {
    assert(m.data != NULL);
    printf("(%d, %d)\n", m.rows, m.cols);
}

int mat_save(Matrix m, const char *path) {
    assert(m.data != NULL && path != NULL);
    FILE *fp = fopen(path, "wb");

    if (fp == NULL) {
        fprintf(stderr, "mat_save: cannot open '%s'\n", path);
        return 0;
    }

    int ele_size = (int) sizeof(MatScalar);
    size_t data_len = (size_t)m.rows * m.cols;

    int is_completed = fwrite(&m.rows, sizeof(int), 1, fp) == 1
                    && fwrite(&m.cols, sizeof(int), 1, fp) == 1
                    && fwrite(&ele_size, sizeof(int), 1, fp) == 1
                    && fwrite(m.data, sizeof(MatScalar), data_len, fp) == data_len;

    if (!is_completed)
        fprintf(stderr, "mat_save: write failed for '%s'\n", path);

    fclose(fp);
    return is_completed;
}

int mat_load(Matrix out, const char *path) {
    assert(out.data != NULL && path != NULL);
    FILE *fp = fopen(path, "rb");

    if (fp == NULL) {
        fprintf(stderr, "mat_load: cannot open '%s'\n", path);
        return 0;
    }

    int rows, cols, ele_size;
    int is_completed = fread(&rows, sizeof(int), 1, fp) == 1
                    && fread(&cols, sizeof(int), 1, fp) == 1
                    && fread(&ele_size, sizeof(int), 1, fp) == 1;

    if (!is_completed)
        fprintf(stderr, "mat_load: cannot read header of '%s'\n", path);
    else if (rows != out.rows || cols != out.cols) {
        fprintf(stderr, "mat_load: '%s' is %d x %d, expected %d x %d\n",
                path, rows, cols, out.rows, out.cols);
        is_completed = 0;
    }
    else if (ele_size != (int)sizeof(MatScalar)) {
        fprintf(stderr, "mat_load: '%s' has %d-byte elements, expected %d\n",
                path, ele_size, (int)sizeof(MatScalar));
        is_completed = 0;
    }
    else {
        size_t data_len = (size_t)out.rows * out.cols;
        is_completed = fread(out.data, ele_size, data_len, fp) == data_len;
        if (!is_completed) fprintf(stderr, "mat_load: reading element values failed for '%s'\n", path);
    }

    fclose(fp);
    return is_completed;
}

MatScalar mat_trace(Matrix m) {
    assert(m.data != NULL);
    assert(m.rows == m.cols);

    double trace = 0.0;
    for (int i = 0; i < m.rows; i++)
        trace += MAT_AT(m, i, i);
    return (MatScalar)trace;
}

MatScalar mat_det(Arena *arena, Matrix m) {
    assert(m.data != NULL);
    assert(m.rows == m.cols);

    size_t mark = arena_mark(arena);

    Matrix work = mat_new(arena, m.rows, m.cols);
    mat_copy(work, m);

    double det = 1.0;
    for (int col = 0; col < work.cols; col++) {
        int pivot_row = col;
        // Find the row where the largest absolute value is inside
        // That row will become a pivot row for one column.
        for (int row = col + 1; row < work.rows; row++) {
            if (fabs(MAT_AT(work, row, col)) > fabs(MAT_AT(work, pivot_row, col)))
                pivot_row = row;
        }

        if (MAT_AT(work, pivot_row, col) == 0) {
            det = 0.0;  // That's mean this matrix is singular
            break;
        }

        // We need a pivot row == col (all pivots should be diagnal), for perform a gaussian elimination
        if (pivot_row != col) {
            mat_swap_rows(work, col, pivot_row);
            det = -det; // each row swap flips the sign of the determinant
        }

        for (int row = col + 1; row < work.rows; row++) {
            MatScalar s = -(MAT_AT(work, row, col) / MAT_AT(work, col, col));
            mat_add_scaled(mat_row_view(work, row), s, mat_row_view(work, col));
        }
        det *= MAT_AT(work, col, col);
    }
    arena_rewind(arena, mark);
    return (MatScalar)det;
}

int mat_inverse(Arena *arena, Matrix out, Matrix m) {
    assert(m.data != NULL && out.data != NULL);
    assert(m.rows == m.cols);
    assert(mat_same_shape(m, out));
    assert(out.data != m.data);

    size_t mark = arena_mark(arena);

    Matrix work = mat_new(arena, m.rows, m.cols);
    mat_copy(work, m);
    mat_identity(out);

    int is_completed = gauss_jordan(work, out);

    arena_rewind(arena, mark);
    return is_completed;
}

int mat_solve(Arena *arena, Matrix x, Matrix A, Matrix b) {
    assert(x.data != NULL && A.data != NULL && b.data != NULL);
    assert(A.rows == A.cols);
    assert(A.rows == b.rows);
    assert(mat_same_shape(b, x));
    assert(x.data != A.data && x.data != b.data);

    size_t mark = arena_mark(arena);

    Matrix work = mat_new(arena, A.rows, A.cols);
    mat_copy(work, A);

    mat_copy(x, b);

    int is_completed = gauss_jordan(work, x);

    arena_rewind(arena, mark);
    return is_completed;
}