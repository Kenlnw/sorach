#include "matrix.h"
#include "lib/arena.h"
#include<stdlib.h>
#include <stdio.h>

#define ARENA_SIZE 1024

int main(void)
{
    Arena arena = arena_create(ARENA_SIZE);

    Matrix eye = mat_new(&arena, 3, 3);
    mat_identity(eye);
    mat_print_shape(eye);
    mat_print(eye);

    Scalar eye_det = mat_det(&arena, eye);
    printf("The determinant of eye is %f\n", eye_det);

    Matrix m = mat_new(&arena, 3, 3);
    mat_copy(m, eye);
    mat_swap_rows(m, 0, 2);
    mat_print_shape(m);
    mat_print(m);

    Scalar m_det = mat_det(&arena, m);
    printf("The determinant of m is %f\n", m_det);

    arena_destroy(&arena);
    return 0;
}
