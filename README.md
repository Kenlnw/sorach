# Sormat

Sora's matrix library: a small matrix library in C99, built on an arena allocator.

- Row-major `Matrix` of `float` (or `double` with `-DMAT_USE_DOUBLE`)
- No hidden allocation: only `mat_new` takes memory from an arena, and every
  other function writes its result into a matrix you pass in
- Elementwise math, broadcasting, reductions along an axis, argmax/argmin
- Matrix product, transpose, determinant, inverse and linear solve
  (Gaussian elimination with partial pivoting)
- Zero-copy row views and reshapes
- Seeded uniform and normal random fill
- Save and load matrices as binary files

## Quick start

```c
#include <stdio.h>
#include <sormat/matrix.h>

int main(void) {
    Arena arena = arena_create(1024 * 1024);

    Matrix a = mat_new(&arena, 2, 3);
    Matrix b = mat_new(&arena, 3, 2);
    mat_from_array(a, (MatScalar[]){1, 2, 3,
                                    4, 5, 6});
    mat_from_array(b, (MatScalar[]){ 7,  8,
                                     9, 10,
                                    11, 12});

    Matrix c = mat_new(&arena, a.rows, b.cols);
    mat_mul(c, a, b);   // c = a @ b
    mat_print(c);

    arena_destroy(&arena);   // frees every matrix at once
    return 0;
}
```

A longer tour is in [`examples/example.c`](examples/example.c).

## Building

You need `gcc` and `make` (on Windows, MinGW/MSYS2 works).

| Command         | What it does                                        |
| --------------- | --------------------------------------------------- |
| `make`          | Build the static library `build/libsormat.a`        |
| `make test`     | Build and run the tests                             |
| `make examples` | Build every program in `examples/` into `build/`    |
| `make clean`    | Delete everything that was built                    |

Options, added to any command above:

- `DOUBLE=1` builds with `double` instead of `float` (output goes to `build_double/`)
- `SANITIZE=1` adds AddressSanitizer and UBSan (Linux and macOS). Run
  `make clean` first when turning it on or off.

### Using the library in your own program

```sh
make
gcc -std=c99 -Iinclude my_program.c build/libsormat.a -lm -o my_program
```

If you build the library with `DOUBLE=1`, compile your program with
`-DMAT_USE_DOUBLE` too and link `build_double/libsormat.a`. Every file in a
program must use the same scalar type.

## How it works

### Arenas

All matrix memory comes from an `Arena` (`include/sormat/arena.h`): one block
allocated up front and handed out in order. Individual matrices are never
freed. Instead you release memory all at once:

- `arena_reset` frees everything
- `arena_mark` / `arena_rewind` free everything allocated since a mark,
  which is handy for temporary matrices inside a loop
- `arena_destroy` gives the block back to the system

Always pass an arena by pointer.

### Conventions

- `MAT_AT(m, r, c)` reads or writes element `(r, c)`. Indexes start at 0.
- A `Matrix` is passed by value, but it only holds a shape and a pointer,
  so copies refer to the same elements.
- Elementwise functions allow the output to be one of the inputs
  (`mat_add(a, a, b)` is `a += b`). Linear algebra functions such as
  `mat_mul` and `mat_transpose` do not.
- Views (`mat_row_view`, `mat_rows_view`, `mat_reshape`) share memory with
  the original matrix. Writing through a view changes the original.
- Axis functions follow NumPy: `MAT_ROW` collapses the rows (one result per
  column), `MAT_COL` collapses the columns (one result per row).
- Caller mistakes such as mismatched shapes stop the program with an
  `assert`. Only `mat_save`, `mat_load`, `mat_inverse` and `mat_solve`
  report failure through their return value.
- Compare floating-point results with `mat_equal(a, b, epsilon)`, not
  exact equality.

The full API, with a comment on every function, is in
[`include/sormat/matrix.h`](include/sormat/matrix.h).

## Project layout

```
include/sormat/   public headers (matrix.h, arena.h)
src/              library source
tests/            unit tests (make test)
examples/         example programs (make examples)
Makefile
```
