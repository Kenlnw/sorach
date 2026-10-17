// example.c: a short tour of the Sorach: SorMat library.
//
// Build and run from the project's top folder:
//     make examples
//     ./build/example

#include <stdio.h>
#include <sorach/matrix.h>

int main(void) {
    // Every matrix lives in an arena. One megabyte is plenty here.
    SorArena arena = sor_arena_create(1024 * 1024);

    // 1. Create matrices and fill them from plain arrays (row by row).
    SorMat a = sor_mat_new(&arena, 2, 3);
    SorMat b = sor_mat_new(&arena, 3, 2);
    sor_mat_from_array(a, (SorScalar[]){1, 2, 3,
                                    4, 5, 6});
    sor_mat_from_array(b, (SorScalar[]){ 7,  8,
                                     9, 10,
                                     11, 12});

    printf("a =\n");
    sor_mat_print(a);
    printf("b =\n");
    sor_mat_print(b);

    // 2. SorMat product. The caller creates the result matrix first;
    //    library functions never allocate on their own.
    SorMat product = sor_mat_new(&arena, a.rows, b.cols);
    sor_mat_mul(product, a, b);
    printf("a @ b =\n");
    sor_mat_print(product);

    // 3. Transpose.
    SorMat a_t = sor_mat_new(&arena, a.cols, a.rows);
    sor_mat_transpose(a_t, a);
    printf("a transposed =\n");
    sor_mat_print(a_t);

    // 4. Elementwise operations may write into one of their inputs.
    sor_mat_scale(a, 2);          // a *= 2
    sor_mat_add(a, a, a);         // a = a + a
    printf("a after scaling by 2 and adding it to itself =\n");
    sor_mat_print(a);

    // 5. Solve a system of linear equations:
    //        2x +  y = 5
    //         x + 3y = 10
    SorMat coeffs = sor_mat_new(&arena, 2, 2);
    SorMat rhs    = sor_mat_new(&arena, 2, 1);
    SorMat x      = sor_mat_new(&arena, 2, 1);
    sor_mat_from_array(coeffs, (SorScalar[]){2, 1,
                                         1, 3});
    sor_mat_from_array(rhs, (SorScalar[]){5, 10});

    // Functions that can fail report it through their return value.
    if (sor_mat_solve(&arena, x, coeffs, rhs)) {
        printf("solution of the system =\n");
        sor_mat_print(x);
    } else {
        printf("the system has no unique solution\n");
    }

    sor_mat_save(x, "x.bin");

    // 6. Free everything at once.
    sor_arena_destroy(&arena);
    return 0;
}