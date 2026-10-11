#ifndef SORMAT_TEST_H
#define SORMAT_TEST_H

#include <stdio.h>

extern int tests_run;
extern int tests_failed;

// Like assert, but records the failure and keeps going.
#define CHECK(cond) do {                                          \
        tests_run++;                                              \
        if (!(cond)) {                                            \
            tests_failed++;                                       \
            printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                         \
    } while (0)

#endif