#include "test.h"

int tests_run = 0;
int tests_failed = 0;

void test_matrix(void);

int main(void) {
    test_matrix();

    printf("%d checks, %d failed\n", tests_run, tests_failed);
    return tests_failed != 0;
}
