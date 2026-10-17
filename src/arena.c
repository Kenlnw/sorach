#include <sorach/arena.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include "check.h"

SorArena sor_arena_create(size_t buff_size) {
    SorArena arena = {
        .buff = malloc(buff_size),
        .buff_size = buff_size,
        .cur_offset = 0
    };

    SOR_CHECK(arena.buff != NULL);
    return arena;
}

void *sor_arena_alloc(SorArena *arena, size_t size, int alignment) {
    SOR_CHECK(arena != NULL && arena->buff != NULL);
    // A power of two has exactly one bit set, so x & (x - 1) is 0.
    SOR_CHECK(alignment > 0 && (alignment & (alignment - 1)) == 0);

    // Align the address itself, not just the offset: malloc only guarantees
    // its own alignment for the start of the buffer (often 16 bytes, 8 on
    // 32-bit Windows), so aligning the offset alone would not give 32 or 64.
    size_t    mask = (size_t)alignment - 1;
    uintptr_t next = (uintptr_t)(arena->buff + arena->cur_offset);
    size_t    pad  = (size_t)((alignment - (next & mask)) & mask);

    // Compare by subtracting, so a huge size cannot wrap around and pass.
    size_t free_bytes = arena->buff_size - arena->cur_offset;
    SOR_CHECK(pad <= free_bytes && size <= free_bytes - pad);   // arena too small

    void *mem = arena->buff + arena->cur_offset + pad;
    arena->cur_offset += pad + size;
    return mem;
}

void *sor_arena_alloc_zero(SorArena *arena, size_t size, int alignment) {
    void *mem = sor_arena_alloc(arena, size, alignment);
    memset(mem, 0, size);
    return mem;
}

void sor_arena_reset(SorArena *arena) {
    SOR_CHECK(arena != NULL);
    arena->cur_offset = 0;
}

size_t sor_arena_mark(SorArena *arena) {
    SOR_CHECK(arena != NULL);
    return arena->cur_offset;
}

void sor_arena_rewind(SorArena *arena, size_t mark) {
    SOR_CHECK(arena != NULL);
    SOR_CHECK(mark <= arena->cur_offset);
    arena->cur_offset = mark;
}

void sor_arena_destroy(SorArena *arena) {
    SOR_CHECK(arena != NULL);
    free(arena->buff);
    arena->buff = NULL;
    arena->buff_size = 0;
    arena->cur_offset = 0;
}