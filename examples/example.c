// example.c: a short tour of the Sormat matrix library.
//
// Build and run from the project's top folder:
//     make examples
//     ./build/example

#include <stdio.h>
#include <sormat/matrix.h>

int main(void) {
    // Every matrix lives in an arena. One megabyte is plenty here.
    Arena arena = arena_create(1024 * 1024);

    // 1. Create matrices and fill them from plain arrays (row by row).
    Matrix a = mat_new(&arena, 2, 3);
    Matrix b = mat_new(&arena, 3, 2);
    mat_from_array(a, (MatScalar[]){1, 2, 3,
                                    4, 5, 6});
    mat_from_array(b, (MatScalar[]){ 7,  8,
                                     9, 10,
                                     11, 12});

    printf("a =\n");
    mat_print(a);
    printf("b =\n");
    mat_print(b);

    // 2. Matrix product. The caller creates the result matrix first;
    //    library functions never allocate on their own.
    Matrix product = mat_new(&arena, a.rows, b.cols);
    mat_mul(product, a, b);
    printf("a @ b =\n");
    mat_print(product);

    // 3. Transpose.
    Matrix a_t = mat_new(&arena, a.cols, a.rows);
    mat_transpose(a_t, a);
    printf("a transposed =\n");
    mat_print(a_t);

    // 4. Elementwise operations may write into one of their inputs.
    mat_scale(a, 2);          // a *= 2
    mat_add(a, a, a);         // a = a + a
    printf("a after scaling by 2 and adding it to itself =\n");
    mat_print(a);

    // 5. Solve a system of linear equations:
    //        2x +  y = 5
    //         x + 3y = 10
    Matrix coeffs = mat_new(&arena, 2, 2);
    Matrix rhs    = mat_new(&arena, 2, 1);
    Matrix x      = mat_new(&arena, 2, 1);
    mat_from_array(coeffs, (MatScalar[]){2, 1,
                                         1, 3});
    mat_from_array(rhs, (MatScalar[]){5, 10});

    // Functions that can fail report it through their return value.
    if (mat_solve(&arena, x, coeffs, rhs)) {
        printf("solution of the system =\n");
        mat_print(x);
    } else {
        printf("the system has no unique solution\n");
    }

    // 6. Free everything at once.
    arena_destroy(&arena);
    return 0;
}