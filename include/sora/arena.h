#ifndef SORA_ARENA_H
#define SORA_ARENA_H

#include <stddef.h>

/*
 * Arena allocator.
 *
 * An arena is one large block of memory obtained up front. Allocations are
 * taken from it in order by advancing an offset, which makes them very cheap.
 * Individual allocations cannot be freed; memory is released all at once with
 * arena_reset, or back to an earlier point with arena_rewind.
 *
 * Typical use:
 *     Arena arena = arena_create(1024 * 1024);
 *     float *x = arena_alloc(&arena, 100 * sizeof(float), ARENA_DEFAULT_ALIGN);
 *     ...
 *     arena_destroy(&arena);
 *
 * Always pass an Arena by pointer. Copying the struct creates two arenas
 * that share one buffer but track their offsets separately.
 */

typedef struct {
    unsigned char *buff;        // start of the memory block
    size_t         buff_size;   // total size of the block, in bytes
    size_t         cur_offset;  // bytes used so far; the next allocation starts here
} Arena;

// Alignment used when the caller has no special requirement. 16 bytes suits
// every standard C type and SSE vector instructions.
#define ARENA_DEFAULT_ALIGN 16

// Create an arena with a buffer of buff_size bytes.
// Stops the program if the memory cannot be obtained.
Arena arena_create(size_t buff_size);

// Allocate size bytes, starting at an offset that is a multiple of alignment.
// The memory is uninitialized. alignment must be greater than 0.
// Stops the program if the arena does not have enough space left.
void *arena_alloc(Arena *arena, size_t size, int alignment);

// Same as arena_alloc, but the returned memory is filled with zeros.
void *arena_alloc_zero(Arena *arena, size_t size, int alignment);

// Release every allocation at once. The buffer is kept and reused, so all
// pointers previously returned from this arena become invalid.
void arena_reset(Arena *arena);

// Return the current position, to be passed to arena_rewind later.
size_t arena_mark(Arena *arena);

// Release everything allocated since mark was taken. Allocations made before
// the mark are kept. Nested marks must be rewound in reverse order.
void arena_rewind(Arena *arena, size_t mark);

// Free the arena's buffer. The arena must not be used afterward.
void arena_destroy(Arena *arena);

#endif