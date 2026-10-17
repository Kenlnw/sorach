/*
 * sorach/arena.h: the arena allocator, part of the Sorach library.
 * https://github.com/Kenlnw/sorach
 * Released under the MIT License; see LICENSE.
 *
 * A SorArena is one large block of memory obtained up front. Allocations are
 * taken from it in order by advancing an offset, which makes them very cheap.
 * Individual allocations cannot be freed. Instead, memory is released all at
 * once with sor_arena_reset, or back to an earlier point with
 * sor_arena_rewind. Every SorMat in Sorach takes its memory from an arena.
 *
 * Conventions:
 *   - Always pass a SorArena by pointer. Copying the struct creates two
 *     arenas that share one buffer but track their offsets separately, so
 *     they would hand out the same memory twice.
 *   - Caller mistakes (an arena that is too small, a bad alignment, a mark
 *     from the future) stop the program with a message naming the failed
 *     check. These checks stay on in every build, including with -DNDEBUG.
 *
 * Typical use:
 *     SorArena arena = sor_arena_create(1024 * 1024);
 *
 *     float *x = sor_arena_alloc(&arena, 100 * sizeof(float),
 *                                SOR_ARENA_DEFAULT_ALIGN);
 *
 *     size_t mark = sor_arena_mark(&arena);
 *     ...                                 // temporary allocations
 *     sor_arena_rewind(&arena, mark);     // frees only the temporaries
 *
 *     sor_arena_destroy(&arena);
 */

#ifndef SORACH_ARENA_H
#define SORACH_ARENA_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

typedef struct {
    unsigned char *buff;         // start of the memory block
    size_t         buff_size;    // total size of the block, in bytes
    size_t         cur_offset;   // bytes used so far; the next allocation starts here
} SorArena;

// Alignment used when the caller has no special requirement. 16 bytes suits
// every standard C type and SSE vector instructions.
#define SOR_ARENA_DEFAULT_ALIGN 16

// ---------------------------------------------------------------------------
// Creating and destroying
// ---------------------------------------------------------------------------

// Create an arena with a buffer of buff_size bytes. buff_size must be
// greater than 0. Stops the program if the memory cannot be obtained.
SorArena sor_arena_create(size_t buff_size);

// Free the arena's buffer. Every pointer and SorMat taken from the arena
// becomes invalid, and the arena must not be used afterward.
void sor_arena_destroy(SorArena *arena);

// ---------------------------------------------------------------------------
// Allocating
// ---------------------------------------------------------------------------

// Allocate size bytes whose address is a multiple of alignment.
// alignment must be a power of two (1, 2, 4, 8, 16, 32, 64, ...); use
// SOR_ARENA_DEFAULT_ALIGN when there is no special requirement. The memory is
// uninitialized. Stops the program if the arena does not have enough space
// left.
void *sor_arena_alloc(SorArena *arena, size_t size, int alignment);

// Same as sor_arena_alloc, but the returned memory is filled with zeros.
void *sor_arena_alloc_zero(SorArena *arena, size_t size, int alignment);

// ---------------------------------------------------------------------------
// Releasing
// Nothing is given back to the system here: the buffer is kept and reused.
// Pointers and SorMats taken from the released part become invalid.
// ---------------------------------------------------------------------------

// Release every allocation at once.
void sor_arena_reset(SorArena *arena);

// Return the current position, to be passed to sor_arena_rewind later.
size_t sor_arena_mark(SorArena *arena);

// Release everything allocated since mark was taken. Allocations made before
// the mark are kept. When marks are nested, rewind them in reverse order:
// the most recent mark first.
void sor_arena_rewind(SorArena *arena, size_t mark);

#ifdef __cplusplus
}
#endif

#endif