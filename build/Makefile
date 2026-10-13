# Makefile for E6E7STAR - closest-point algorithms for the E6* and E7* lattices
#
# Usage:
#   make              - build the test executable
#   make test         - build and run the test program
#   make clean        - remove build artifacts
#   make rebuild      - clean and build
#   make info         - show the build configuration
#   make help         - show this help message
#
# Variables:
#   CC         - C compiler (default: gcc)
#   CFLAGS     - additional compiler flags, appended to the built-in ones
#   DEBUG      - set to 1 for debug build with -g -O0
#
# Built-in flags: -std=c89 -pedantic -Wall -Wextra -O2, and the options that turn off automatic vectorization:
#   gcc:    -fno-tree-vectorize
#   clang:  -fno-vectorize -fno-slp-vectorize
# No -march option is given. So the program is plain C89, compiled to scalar code for the baseline instruction
# set of the target, and runs in one thread. All timings of the project refer to this configuration.

# ============================================================================
# Configuration
# ============================================================================

# Directories
SRC_DIR     := ../src
INC_DIR     := ../include
BUILD_DIR   := .
BIN_DIR     := $(BUILD_DIR)
OBJ_DIR     := $(BUILD_DIR)/obj

# Source files
SOURCES     := $(SRC_DIR)/e6e7star_test.c \
              $(SRC_DIR)/e6e7star_CS.c \
              $(SRC_DIR)/e6e7star_TYK.c \
              $(SRC_DIR)/e6e7star_CS84.c \
              $(SRC_DIR)/e6e7star_MSC.c \
              $(SRC_DIR)/e6e7star_SE.c \
              $(SRC_DIR)/e6e7star_CA.c \
              $(SRC_DIR)/e6e7star_A2.c \
              $(SRC_DIR)/e6e7star_CS_opt.c \
              $(SRC_DIR)/e6e7star_CS84_opt.c \
              $(SRC_DIR)/e6e7star_sweep.c \
              $(SRC_DIR)/e6e7star_sweep_opt.c \
              $(SRC_DIR)/e6e7star_net.c \
              $(SRC_DIR)/e6e7star_bucket.c \
              $(SRC_DIR)/e6e7star_Rn.c

# Object files
OBJECTS     := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

# Executable
TARGET      := $(BIN_DIR)/e6e7star_test

# Compiler: gcc, unless CC is set on the command line or in the environment
ifeq ($(origin CC),default)
CC          := gcc
endif

# Built-in flags
CFLAGS_BASE := -std=c89 -pedantic -Wall -Wextra
CFLAGS_OPT  := -O2
CFLAGS_DBG  := -g -O0 -DDEBUG

# Scalar code only: turn off the automatic vectorizers of the compiler. These options must follow -O2.
# "make CFLAGS=..." can turn them on again, or add options such as -march=native.
CC_VERSION  := $(shell $(CC) --version 2>/dev/null)
ifneq (,$(findstring clang,$(CC_VERSION)))
CFLAGS_OPT  += -fno-vectorize -fno-slp-vectorize
else ifneq (,$(findstring Free Software Foundation,$(CC_VERSION)))
CFLAGS_OPT  += -fno-tree-vectorize
endif

# Determine optimization flags based on DEBUG variable
ifeq ($(DEBUG),1)
    CFLAGS_MODE := $(CFLAGS_DBG)
else
    CFLAGS_MODE := $(CFLAGS_OPT)
endif

# All compiler flags. CFLAGS given by the user come last, so that they can override the built-in ones.
ALL_CFLAGS  := $(CFLAGS_BASE) $(CFLAGS_MODE) -I$(INC_DIR) $(CFLAGS)

# Linker flags
LDFLAGS     := -lm

# ============================================================================
# Build Rules
# ============================================================================

# Default target
.PHONY: all
all: $(TARGET)

# Create object directory
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/e6e7star.h | $(OBJ_DIR)
	$(CC) $(ALL_CFLAGS) -c $< -o $@

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CC) $(ALL_CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# ============================================================================
# Maintenance Targets
# ============================================================================

.PHONY: clean
clean:
	@rm -rf $(OBJ_DIR)
	@rm -f $(TARGET)
	@echo "Clean complete"

# Clean and build, one after the other (also under make -j)
.PHONY: rebuild
rebuild:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory all

# Build and run the test program with its default parameters
.PHONY: test
test: $(TARGET)
	$(TARGET)

.PHONY: help
help:
	@echo "E6E7STAR Build System"
	@echo "====================="
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build the test executable (default)"
	@echo "  make test         - Build and run the test program"
	@echo "  make clean        - Remove all build artifacts"
	@echo "  make rebuild      - Clean and build"
	@echo "  make info         - Show the build configuration"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  CC=<compiler>     - Specify C compiler (default: gcc)"
	@echo "  CFLAGS=<flags>    - Specify additional compiler flags"
	@echo "  DEBUG=1           - Build with debug symbols and no optimization"
	@echo ""
	@echo "Examples:"
	@echo "  make                          # Standard optimized build"
	@echo "  make DEBUG=1                  # Debug build"
	@echo "  make CC=clang                 # Use clang instead of gcc"
	@echo "  make CFLAGS=-march=native     # Add a compiler flag"
	@echo "  make rebuild                  # Clean and rebuild"
	@echo ""
	@echo "Output:"
	@echo "  Executable: $(TARGET)"
	@echo "  Objects:    $(OBJ_DIR)/"

.PHONY: info
info:
	@echo "Build Configuration:"
	@echo "  Compiler:    $(CC)"
	@echo "  CFLAGS:      $(ALL_CFLAGS)"
	@echo "  LDFLAGS:     $(LDFLAGS)"
	@echo "  Source dir:  $(SRC_DIR)"
	@echo "  Include dir: $(INC_DIR)"
	@echo "  Build dir:   $(BUILD_DIR)"
	@echo "  Target:      $(TARGET)"
