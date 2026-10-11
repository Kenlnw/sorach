# Makefile for the Sormat library
#
#   make            build the static library (build/libsormat.a)
#   make test       build and run the tests
#   make examples   build the programs in examples/ (into the build folder)
#   make clean      delete everything that was built
#
# Options (add to any command above):
#   DOUBLE=1        build with double instead of float (uses build_double/)
#   SANITIZE=1      add AddressSanitizer and UBSan (Linux and macOS);
#                   run "make clean" first when turning this on or off

CC     := gcc
CFLAGS := -std=c99 -Wall -Wextra -g -Iinclude
LDLIBS := -lm
BUILD  := build

ifdef DOUBLE
CFLAGS += -DMAT_USE_DOUBLE
BUILD  := build_double
endif

ifdef SANITIZE
CFLAGS += -fsanitize=address,undefined
endif

ifeq ($(OS),Windows_NT)
EXE := .exe
endif

LIB      := $(BUILD)/libsormat.a
TEST_BIN := $(BUILD)/run_tests$(EXE)

# Every .c file in src/ goes into the library; every .c file in tests/ goes
# into the test program. New files are picked up automatically.
# All object files share one folder, so a file in tests/ must not have the
# same name as a file in src/.
LIB_SRC  := $(wildcard src/*.c)
TEST_SRC := $(wildcard tests/*.c)
LIB_OBJ  := $(patsubst src/%.c,$(BUILD)/%.o,$(LIB_SRC))
TEST_OBJ := $(patsubst tests/%.c,$(BUILD)/%.o,$(TEST_SRC))

# Each .c file in examples/ becomes its own program.
EX_SRC   := $(wildcard examples/*.c)
EX_BIN   := $(patsubst examples/%.c,$(BUILD)/%$(EXE),$(EX_SRC))

.PHONY: all test examples clean

all: $(LIB)

# Pack the library's object files into a static library.
$(LIB): $(LIB_OBJ)
	$(AR) rcs $@ $^

# Link the test program against the library.
$(TEST_BIN): $(TEST_OBJ) $(LIB)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

test: $(TEST_BIN)
	./$(TEST_BIN)

examples: $(EX_BIN)

# Compile one example and link it against the library.
$(BUILD)/%$(EXE): examples/%.c $(LIB) | $(BUILD)
	$(CC) $(CFLAGS) $< $(LIB) -o $@ $(LDLIBS)

# Compile one .c file into one .o file. -MMD -MP also writes a .d file
# listing the headers it includes, so editing a header triggers a rebuild.
$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/%.o: tests/%.c | $(BUILD)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD):
	mkdir $(BUILD)

clean:
	rm -rf build build_double

-include $(wildcard $(BUILD)/*.d)