#ifndef SORACH_CHECK_H
#define SORACH_CHECK_H

// Internal to the library: not installed with the public headers.

#include <stdio.h>
#include <stdlib.h>

// Like assert, but never compiled out: -DNDEBUG does not remove it.
// If cond is false, print where it failed and stop the program.
#define SOR_CHECK(cond)                                                   \
    do {                                                                  \
        if (!(cond)) {                                                    \
            fprintf(stderr, "sorach: %s:%d: %s: check failed: %s\n",      \
                    __FILE__, __LINE__, __func__, #cond);                 \
            abort();                                                      \
        }                                                                 \
    } while (0)

#endif