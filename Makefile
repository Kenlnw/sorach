# Makefile for the Sorach library
#
#   make            build the static library (build/libsorach.a)
#   make test       build and run the tests
#   make examples   build the programs in examples/ (into the build folder)
#   make install    copy the headers and the library under PREFIX
#   make uninstall  remove what make install copied
#                   (install and uninstall need a Unix-style shell: Linux,
#                   macOS, or the MSYS2 / Git Bash terminal on Windows)
#   make clean      delete everything that was built
#
# Options (add to any command above):
#   DOUBLE=1        use double instead of float
#   RELEASE=1       optimize (-O2) instead of building for debugging (-g);
#                   the library's checks stay on either way
#   SANITIZE=1      add AddressSanitizer and UBSan (on Windows: UBSan only,
#                   in trap mode; see below)
#   PREFIX=path     where make install puts files (default /usr/local)
#
# Each combination of DOUBLE, RELEASE and SANITIZE builds into its own
# folder (build, build_double, build_release, build_double_sanitize, ...),
# so switching options never mixes files built with different settings.

CC     := gcc
CFLAGS := -std=c99 -Wall -Wextra -Iinclude
LDLIBS := -lm
BUILD  := build

ifdef DOUBLE
CFLAGS += -DSOR_USE_DOUBLE
BUILD  := $(BUILD)_double
endif

ifdef RELEASE
CFLAGS += -O2
BUILD  := $(BUILD)_release
else
CFLAGS += -g
endif

# MinGW gcc on Windows has no AddressSanitizer or UBSan runtime libraries,
# so there SANITIZE=1 uses UBSan's trap mode, which needs no library: the
# program stops at the first undefined operation instead of printing a report.
# It does not detect memory errors; for those, build on Linux, macOS or WSL.
ifdef SANITIZE
ifeq ($(OS),Windows_NT)
CFLAGS += -fsanitize=undefined -fsanitize-undefined-trap-on-error
else
CFLAGS += -fsanitize=address,undefined
endif
BUILD  := $(BUILD)_sanitize
endif

ifeq ($(OS),Windows_NT)
EXE := .exe
endif

# How to delete a folder. When make runs its commands through Windows' cmd
# (for example mingw32-make started from PowerShell) there is no rm, so use
# rmdir instead. cmd's echo keeps quotes and a Unix shell's does not, which
# tells the two apart.
ifeq ($(shell echo "x"),"x")
RMDIR := rmdir /s /q
else
RMDIR := rm -rf
endif

LIB      := $(BUILD)/libsorach.a
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

# Where make install copies files.
PREFIX     ?= /usr/local
INCLUDEDIR := $(PREFIX)/include
LIBDIR     := $(PREFIX)/lib

.PHONY: all test examples install uninstall clean

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

# Only the public headers in include/sorach/ are installed; src/check.h is
# internal to the library and stays behind.
install: $(LIB)
	mkdir -p $(DESTDIR)$(INCLUDEDIR)/sorach $(DESTDIR)$(LIBDIR)
	cp include/sorach/*.h $(DESTDIR)$(INCLUDEDIR)/sorach/
	cp $(LIB) $(DESTDIR)$(LIBDIR)/

uninstall:
	rm -rf $(DESTDIR)$(INCLUDEDIR)/sorach
	rm -f $(DESTDIR)$(LIBDIR)/libsorach.a

# Deletes build/ and every build_* folder that exists.
clean:
	$(foreach d,$(wildcard build build_*),$(RMDIR) $(d) &&) echo cleaned

-include $(wildcard $(BUILD)/*.d)