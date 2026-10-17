#ifndef SORACH_CORE_H
#define SORACH_CORE_H

#define SORACH_VERSION_MAJOR 0
#define SORACH_VERSION_MINOR 1
#define SORACH_VERSION_PATCH 0

// The version as one number, for comparisons in the preprocessor:
//     #if SORACH_VERSION >= SORACH_MAKE_VERSION(0, 2, 0)
#define SORACH_MAKE_VERSION(major, minor, patch) \
    ((major) * 10000 + (minor) * 100 + (patch))
#define SORACH_VERSION \
    SORACH_MAKE_VERSION(SORACH_VERSION_MAJOR, SORACH_VERSION_MINOR, SORACH_VERSION_PATCH)

// The element type used by every part of Sorach. float by default; compile
// with -DSORACH_USE_DOUBLE to use double instead. Every file in a program,
// including the library itself, must use the same setting.
#ifdef SORACH_USE_DOUBLE
typedef double SorScalar;
#else
typedef float SorScalar;
#endif

#endif