/*
 * sorach/sorach.h: includes every part of the Sorach library.
 * https://github.com/Kenlnw/sorach
 * Released under the MIT License; see LICENSE.
 *
 * Most programs only need this one line:
 *     #include <sorach/sorach.h>
 *
 * Each part can also be included on its own, for example
 * <sorach/matrix.h>, if a program uses only some of the library.
 */

#ifndef SORACH_H
#define SORACH_H

#include "core.h"     // SorScalar, version, build settings
#include "arena.h"    // SorArena: the memory every other part allocates from
#include "matrix.h"   // SorMat: 2-D matrices

#endif