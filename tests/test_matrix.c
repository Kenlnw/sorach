#include <sorach/matrix.h>
#include "test.h"

#define EPS ((SorScalar)1e-5)   //Epsilon

// Make a matrix and fill it from an array in one step.
static SorMat make(SorArena *arena, int rows, int cols, const SorScalar *values) {
    SorMat m = sor_mat_new(arena, rows, cols);
    sor_mat_from_array(m, values);
    return m;
}

static void test_mul(SorArena *arena) {
    SorMat a = make(arena, 2, 2, (SorScalar[]){1, 2, 3, 4});
    SorMat b = make(arena, 2, 2, (SorScalar[]){5, 6, 7, 8});
    SorMat want = make(arena, 2, 2, (SorScalar[]){19, 22, 43, 50});
    SorMat got = sor_mat_new(arena, 2, 2);

    sor_mat_mul(got, a, b);
    CHECK(sor_mat_equal(got, want, EPS));
}

static void test_same_shape(SorArena *arena) {
    SorMat a = sor_mat_new(arena, 2, 3);
    CHECK(sor_mat_same_shape(a, sor_mat_new(arena, 2, 3)) == 1);
    CHECK(sor_mat_same_shape(a, sor_mat_new(arena, 3, 3)) == 0);   // rows differ
    CHECK(sor_mat_same_shape(a, sor_mat_new(arena, 2, 5)) == 0);   // cols differ
}

static void test_add_vec_column(SorArena *arena) {
    SorMat m = make(arena, 2, 3, (SorScalar[]){1, 2, 3, 4, 5, 6});
    SorMat v = make(arena, 2, 1, (SorScalar[]){10, 20});
    SorMat want = make(arena, 2, 3, (SorScalar[]){11, 12, 13, 24, 25, 26});

    sor_mat_add_vec(m, v);
    CHECK(sor_mat_equal(m, want, EPS));
}

void test_matrix(void) {
    SorArena arena = sor_arena_create(1024 * 1024);

    test_mul(&arena);
    sor_arena_reset(&arena);

    test_same_shape(&arena);
    sor_arena_reset(&arena);

    test_add_vec_column(&arena);
    sor_arena_reset(&arena);

    sor_arena_destroy(&arena);
}