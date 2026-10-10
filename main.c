#include "matrix.h"
#include "lib/arena.h"
#include<stdlib.h>
#include <stdio.h>

#define ARENA_SIZE 1024

int main(void)
{
    Arena arena = arena_create(ARENA_SIZE);

    Matrix c = mat_new(&arena, 2, 2);
    const char *src = "c.bin";
    mat_load(c, src);
    mat_print(c);
    printf("\n");

    int *c_argmin = arena_alloc(&arena, c.rows, DEFAULT_ALIGN);
    mat_argmin(c_argmin, c.rows, c, ROW);

    for (int i = 0; i < c.rows; i++)
        printf("%d ", c_argmin[i]);
    printf("\n");

    arena_destroy(&arena);
    return 0;
}
