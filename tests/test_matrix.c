#include <sora/matrix.h>
#include "test.h"

#define EPS ((MatScalar)1e-5)   //Epsilon

// Make a matrix and fill it from an array in one step.
static Matrix make(Arena *arena, int rows, int cols, const MatScalar *values) {
    Matrix m = mat_new(arena, rows, cols);
    mat_from_array(m, values);
    return m;
}

static void test_mul(Arena *arena) {
    Matrix a = make(arena, 2, 2, (MatScalar[]){1, 2, 3, 4});
    Matrix b = make(arena, 2, 2, (MatScalar[]){5, 6, 7, 8});
    Matrix want = make(arena, 2, 2, (MatScalar[]){19, 22, 43, 50});
    Matrix got = mat_new(arena, 2, 2);

    mat_mul(got, a, b);
    CHECK(mat_equal(got, want, EPS));
}

static void test_same_shape(Arena *arena) {
    Matrix a = mat_new(arena, 2, 3);
    CHECK(mat_same_shape(a, mat_new(arena, 2, 3)) == 1);
    CHECK(mat_same_shape(a, mat_new(arena, 3, 3)) == 0);   // rows differ
    CHECK(mat_same_shape(a, mat_new(arena, 2, 5)) == 0);   // cols differ
}

static void test_add_vec_column(Arena *arena) {
    Matrix m = make(arena, 2, 3, (MatScalar[]){1, 2, 3, 4, 5, 6});
    Matrix v = make(arena, 2, 1, (MatScalar[]){10, 20});
    Matrix want = make(arena, 2, 3, (MatScalar[]){11, 12, 13, 24, 25, 26});

    mat_add_vec(m, v);
    CHECK(mat_equal(m, want, EPS));
}

void test_matrix(void) {
    Arena arena = arena_create(1024 * 1024);

    test_mul(&arena);
    arena_reset(&arena);

    test_same_shape(&arena);
    arena_reset(&arena);

    test_add_vec_column(&arena);
    arena_reset(&arena);

    arena_destroy(&arena);
}