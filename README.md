# Sorach

Sorach is a small linear algebra library for C99. It gives you dense
matrices of `float` or `double`, the operations you need to work with them,
and an arena allocator that keeps memory management simple and predictable.

- Row-major `SorMat` matrices of `float`, or `double` with `-DSOR_USE_DOUBLE`
- No hidden allocation: only `sor_mat_new` takes memory, from an arena you
  create, and every other function writes its result into a matrix you pass in
- Matrix product (also with the first matrix transposed), transpose, dot
  product, trace, diagonal and identity matrices
- Determinant, inverse and linear solve for square systems, by Gaussian
  elimination with partial pivoting
- Elementwise math, broadcasting a vector across rows or columns, and
  reductions along an axis: sum, mean, standard deviation, max, min,
  argmax and argmin
- Zero-copy row views and reshapes
- Seeded uniform and normal random fill
- Save and load matrices as binary files, one per file or several in one file
- Shape and index checks that stay on in release builds

## Quick start

```c
#include <stdio.h>
#include <sorach/sorach.h>

int main(void) {
    SorArena arena = sor_arena_create(1024 * 1024);

    SorMat a = sor_mat_new(&arena, 2, 3);
    SorMat b = sor_mat_new(&arena, 3, 2);
    sor_mat_from_array(a, (SorScalar[]){1, 2, 3,
                                        4, 5, 6});
    sor_mat_from_array(b, (SorScalar[]){ 7,  8,
                                         9, 10,
                                        11, 12});

    SorMat c = sor_mat_new(&arena, a.rows, b.cols);
    sor_mat_mul(c, a, b);   // c = a @ b
    sor_mat_print(c);

    sor_arena_destroy(&arena);   // frees every matrix at once
    return 0;
}
```

Solving a system of linear equations:

```c
//  2x +  y = 5
//   x + 3y = 10
SorMat A = sor_mat_new(&arena, 2, 2);
SorMat b = sor_mat_new(&arena, 2, 1);
SorMat x = sor_mat_new(&arena, 2, 1);
sor_mat_from_array(A, (SorScalar[]){2, 1,
                                    1, 3});
sor_mat_from_array(b, (SorScalar[]){5, 10});

if (sor_mat_solve(&arena, x, A, b))
    sor_mat_print(x);   // x = 1, y = 3
else
    printf("the system has no unique solution\n");
```

A longer tour is in [`examples/example.c`](examples/example.c).

## Building

You need `gcc` and `make`. On Windows, use MinGW through MSYS2.

| Command          | What it does                                             |
| ---------------- | -------------------------------------------------------- |
| `make`           | Build the static library `build/libsorach.a`             |
| `make test`      | Build and run the tests                                  |
| `make examples`  | Build every program in `examples/` into the build folder |
| `make install`   | Copy the headers and the library under `PREFIX`          |
| `make uninstall` | Remove what `make install` copied                        |
| `make clean`     | Delete everything that was built                         |

Options, added to any command above:

| Option        | Effect                                                          |
| ------------- | --------------------------------------------------------------- |
| `DOUBLE=1`    | Use `double` instead of `float`                                 |
| `RELEASE=1`   | Optimize with `-O2` instead of building for debugging           |
| `SANITIZE=1`  | Add AddressSanitizer and UBSan (Linux and macOS)                |
| `PREFIX=path` | Where `make install` puts files (default `/usr/local`)          |

Each combination of `DOUBLE`, `RELEASE` and `SANITIZE` builds into its own
folder, such as `build/`, `build_double/` or `build_double_release/`, so
switching options never mixes files built with different settings.

### Using the library in your own program

Either install it:

```sh
make RELEASE=1
sudo make install RELEASE=1
gcc -std=c99 my_program.c -lsorach -lm -o my_program
```

or link it straight from the build folder without installing:

```sh
make
gcc -std=c99 -Iinclude my_program.c build/libsorach.a -lm -o my_program
```

If you build the library with `DOUBLE=1`, compile your own program with
`-DSOR_USE_DOUBLE` too, and link the library from the matching build folder.
Every file in a program must use the same scalar type.

## How it works

### Arenas

All matrix memory comes from a `SorArena` (`include/sorach/arena.h`): one
block allocated up front and handed out in order. Individual matrices are
never freed. Instead you release memory all at once:

- `sor_arena_reset` frees everything
- `sor_arena_mark` / `sor_arena_rewind` free everything allocated since a
  mark, which is handy for temporary matrices inside a loop
- `sor_arena_destroy` gives the block back to the system

```c
for (int step = 0; step < steps; step++) {
    size_t mark = sor_arena_mark(&arena);

    SorMat tmp = sor_mat_new(&arena, n, n);   // only needed for this step
    ...

    sor_arena_rewind(&arena, mark);           // tmp is freed here
}
```

Always pass an arena by pointer.

### Conventions

- `SOR_MAT_AT(m, r, c)` reads or writes element `(r, c)`. Indexes start at 0.
- A `SorMat` is passed by value, but it only holds a shape and a pointer,
  so copies refer to the same elements.
- Elementwise functions allow the output to be one of the inputs
  (`sor_mat_add(a, a, b)` is `a += b`). Linear algebra functions such as
  `sor_mat_mul` and `sor_mat_transpose` do not.
- Views (`sor_mat_row_view`, `sor_mat_rows_view`, `sor_mat_reshape`) share
  memory with the original matrix. Writing through a view changes the
  original.
- Axis functions follow NumPy, and the axis values are NumPy's axis numbers:
  `SOR_MAT_ROW` (axis 0) collapses the rows, giving one result per column, and
  `SOR_MAT_COL` (axis 1) collapses the columns, giving one result per row.
- Compare floating-point results with `sor_mat_equal(a, b, epsilon)`, not
  exact equality.

### Errors

Caller mistakes, such as mismatched shapes or an index out of range, stop
the program with a message naming the file, line and failed check:

```
sorach: src/matrix.c:222: sor_mat_mul: check failed: a.cols == b.rows
```

These checks stay on in every build, including with `-DNDEBUG`, so a mistake
never turns into silent memory corruption.

Failures that can happen in a correct program are reported through the
return value instead: 1 for success and 0 for failure. These are the file
functions and the functions that can meet a singular matrix
(`sor_mat_inverse` and `sor_mat_solve`).

### Saving and loading

`sor_mat_save` and `sor_mat_load` store one matrix per file. To store several
matrices in one file, such as all the weights of a model, open the file
yourself and use the stream functions:

```c
FILE *fp = fopen("model.bin", "wb");
sor_mat_write(w1, fp);
sor_mat_write(w2, fp);
fclose(fp);

fp = fopen("model.bin", "rb");
sor_mat_read_into(w1, fp);   // into a matrix that already exists
sor_mat_read(&arena, &w3, fp);   // or create a new one in the arena
fclose(fp);
```

Always open the file in binary mode (`"wb"` or `"rb"`).

Each stored matrix starts with a short header: the bytes `SMAT`, a format
version, the shape and the element size. Loading checks all of these, so a
file of the wrong kind, the wrong shape, or the wrong scalar type is
reported instead of read as garbage. Numbers are stored in the byte order of
the machine that wrote the file, so a file cannot be moved between
little-endian and big-endian machines; this is detected and reported too.

## API reference

The full API, with a comment on every function, is in the headers:

| Header                                                 | Contents                                  |
| ------------------------------------------------------ | ----------------------------------------- |
| [`sorach/sorach.h`](include/sorach/sorach.h)           | Includes everything below                 |
| [`sorach/core.h`](include/sorach/core.h)               | `SorScalar`, version, build settings      |
| [`sorach/arena.h`](include/sorach/arena.h)             | `SorArena`, the memory allocator          |
| [`sorach/matrix.h`](include/sorach/matrix.h)           | `SorMat` and every matrix function        |

## Project layout

```
include/sorach/   public headers
src/              library source
tests/            unit tests (make test)
examples/         example programs (make examples)
LICENSE
Makefile
```

## License

Sorach is released under the MIT License. See [`LICENSE`](LICENSE).