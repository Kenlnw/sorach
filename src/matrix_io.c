#include <sorach/matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include "check.h"

static const char     SOR_MAT_FILE_MAGIC[4] = {'S', 'M', 'A', 'T'};
static const uint32_t SOR_MAT_FILE_VERSION  = 1;

// Header: magic (4 bytes), version, rows, cols, element size (4 bytes each).
static int write_header(SorMat m, FILE *fp) {
    int32_t rows = m.rows, cols = m.cols, elem = (int32_t)sizeof(SorScalar);
    return fwrite(SOR_MAT_FILE_MAGIC, 1, 4, fp) == 4
        && fwrite(&SOR_MAT_FILE_VERSION, sizeof SOR_MAT_FILE_VERSION, 1, fp) == 1
        && fwrite(&rows, sizeof rows, 1, fp) == 1
        && fwrite(&cols, sizeof cols, 1, fp) == 1
        && fwrite(&elem, sizeof elem, 1, fp) == 1;
}

// Read a header and check every field. On success, sets *rows and *cols to
// a valid shape and returns 1. On failure, prints why and returns 0.
// who is the public function's name, used in the messages.
static int read_header(FILE *fp, int32_t *rows, int32_t *cols, const char *who) {
    char     magic[4];
    uint32_t version;
    int32_t  elem;

    if (fread(magic, 1, 4, fp) != 4
        || fread(&version, sizeof version, 1, fp) != 1
        || fread(rows, sizeof *rows, 1, fp) != 1
        || fread(cols, sizeof *cols, 1, fp) != 1
        || fread(&elem, sizeof elem, 1, fp) != 1) {
        fprintf(stderr, "%s: cannot read the header\n", who);
        return 0;
    }
    if (memcmp(magic, SOR_MAT_FILE_MAGIC, 4) != 0) {
        fprintf(stderr, "%s: not a sorach matrix (wrong magic bytes)\n", who);
        return 0;
    }
    if (version != SOR_MAT_FILE_VERSION) {
        fprintf(stderr, "%s: format version %lu is not supported (expected %lu)\n",
                who, (unsigned long)version, (unsigned long)SOR_MAT_FILE_VERSION);
        return 0;
    }
    // Checked before any element is read, so fread below always uses
    // sizeof(SorScalar) and never a size taken from the file.
    if (elem != (int32_t)sizeof(SorScalar)) {
        fprintf(stderr, "%s: file has %ld-byte elements, this build uses %d\n",
                who, (long)elem, (int)sizeof(SorScalar));
        return 0;
    }
    // A damaged file can hold any numbers here. Refuse shapes that
    // sor_mat_new would reject, or whose size does not fit in an int.
    if (*rows <= 0 || *cols <= 0 || *rows > INT_MAX / *cols) {
        fprintf(stderr, "%s: invalid shape %ld x %ld\n", who, (long)*rows, (long)*cols);
        return 0;
    }
    return 1;
}

static int read_data(SorMat out, FILE *fp, const char *who) {
    size_t n = (size_t)out.rows * out.cols;
    if (fread(out.data, sizeof(SorScalar), n, fp) != n) {
        fprintf(stderr, "%s: the file ends before all elements were read\n", who);
        return 0;
    }
    return 1;
}

int sor_mat_write(SorMat m, FILE *fp) {
    SOR_CHECK(m.data != NULL && fp != NULL);
    size_t n = (size_t)m.rows * m.cols;
    if (!write_header(m, fp) || fwrite(m.data, sizeof(SorScalar), n, fp) != n) {
        fprintf(stderr, "sor_mat_write: write failed\n");
        return 0;
    }
    return 1;
}

int sor_mat_read_into(SorMat out, FILE *fp) {
    SOR_CHECK(out.data != NULL && fp != NULL);
    int32_t rows, cols;
    if (!read_header(fp, &rows, &cols, "sor_mat_read_into"))
        return 0;
    if (rows != out.rows || cols != out.cols) {
        fprintf(stderr, "sor_mat_read_into: stored matrix is %ld x %ld, expected %d x %d\n",
                (long)rows, (long)cols, out.rows, out.cols);
        return 0;
    }
    return read_data(out, fp, "sor_mat_read_into");
}

int sor_mat_read(SorArena *arena, SorMat *out, FILE *fp) {
    SOR_CHECK(arena != NULL && out != NULL && fp != NULL);
    int32_t rows, cols;
    if (!read_header(fp, &rows, &cols, "sor_mat_read"))
        return 0;

    // Check the size before allocating: a damaged file claiming a huge
    // matrix must be refused here, not stop the program inside sor_mat_new.
    size_t need_size = (size_t)rows * cols * sizeof(SorScalar) + SOR_ARENA_DEFAULT_ALIGN;
    if (need_size > arena->buff_size - arena->cur_offset) {
        fprintf(stderr, "sor_mat_read: a %ld x %ld matrix does not fit in the arena\n",
                (long)rows, (long)cols);
        return 0;
    }

    size_t mark = sor_arena_mark(arena);
    SorMat m = sor_mat_new(arena, rows, cols);
    if (!read_data(m, fp, "sor_mat_read")) {
        sor_arena_rewind(arena, mark);
        return 0;
    }
    *out = m;
    return 1;
}

int sor_mat_save(SorMat m, const char *path) {
    SOR_CHECK(m.data != NULL && path != NULL);
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        fprintf(stderr, "sor_mat_save: cannot open '%s'\n", path);
        return 0;
    }
    int is_completed = sor_mat_write(m, fp);
    // fwrite sends the bytes in a buffer insize FILE. In some case, the fclose will finish before fwrite.
    // The last bytes may not have reached the file so we need to error out the unfinished fwrite
    if (fclose(fp) != 0) is_completed = 0;
    if (!is_completed) fprintf(stderr, "sor_mat_save: failed to save '%s'\n", path);
    return is_completed;
}

int sor_mat_load(SorMat out, const char *path) {
    SOR_CHECK(out.data != NULL && path != NULL);
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) {
        fprintf(stderr, "sor_mat_load: cannot open '%s'\n", path);
        return 0;
    }
    int ok = sor_mat_read_into(out, fp);
    fclose(fp);   // fread doesn't use the buffer so fread will end immediately. fclose will always finish after fread
    if (!ok) fprintf(stderr, "sor_mat_load: failed to load '%s'\n", path);
    return ok;
}