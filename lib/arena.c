#include "arena.h"
#include<string.h>
#include<stdlib.h>
#include<assert.h>

Arena arena_create(size_t buff_size) {
    Arena arena = {
        .buff = malloc(buff_size),
        .buff_size = buff_size,
        .cur_offset = 0
    };

    assert(arena.buff != NULL);
    return arena;
}

void *arena_alloc(Arena *arena, size_t size, int alignment) {
    assert(alignment > 0);  // check that alignment > 0 or not. If not, stop program
    size_t offset = arena->cur_offset;
    size_t rem = offset % alignment;

    if (rem > 0) offset += alignment - rem;

    // check that the allocate memory will be not overflow
    assert(offset + size <= arena->buff_size);

    // allocate the memory in the buffer (*mem at the starting point of allocated memory)
    void *mem = arena->buff + offset;
    arena->cur_offset = offset + size;

    return mem;
}

void *arena_alloc_zero(Arena *arena, size_t size, int alignment) {
    void *mem = arena_alloc(arena, size, alignment);
    memset(mem, 0, size);
    return mem;
}

void arena_reset(Arena *arena) {
    arena->cur_offset = 0;
}

size_t arena_mark(Arena *arena) {
    return arena->cur_offset;
}

void arena_rewind(Arena *arena, size_t mark) {
    assert(mark <= arena->cur_offset);
    arena->cur_offset = mark;
}

void arena_destroy(Arena *arena) {
    free(arena->buff);
    arena->buff = NULL;
    arena->buff_size = 0;
    arena->cur_offset = 0;
}