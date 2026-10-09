#ifndef ARENA_H
#define ARENA_H

#include<stddef.h>

typedef struct {
    unsigned char *buff;
    size_t buff_size;
    size_t cur_offset;
} Arena;

#define DEFAULT_ALIGN 16

Arena arena_create(size_t buff_size);
void *arena_alloc(Arena *arena, size_t size, int alignment);
void *arena_alloc_zero(Arena *arena, size_t size, int alignment);
void arena_reset(Arena *arena);
void arena_destroy(Arena *arena);

#endif