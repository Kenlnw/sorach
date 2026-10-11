// example_1.c: Create matrices and fill them from plain arrays (row by row) then print them.

#include <stdio.h>
#include <sora/matrix.h>

int main(void) {
    // Every matrix lives in an arena. One megabyte is plenty here.
    Arena arena = arena_create(1024 * 1024);

    Matrix a = mat_new(&arena, 2, 3);
    mat_from_array(a, (MatScalar[]){1, 2, 3,
                                    4, 5, 6});

    Matrix b = mat_new(&arena, 3, 2);
    mat_from_array(b, (MatScalar[]){7, 8,
                                    9, 10,
                                    11, 12});

   printf("a =\n");
   mat_print(a);

   printf("b =\n");
   mat_print(b);

    return 0;
}