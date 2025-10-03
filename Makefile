# Compiler
CC = gcc
MPICC = mpicc

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# Include paths
INCLUDES = -I$(INCLUDE_DIR)

# Flags
MATH = -lm
CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -Wshadow -Wuninitialized -W -Wno-sign-compare -Wno-pedantic -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-unused-variable
OPT_FLAGS = -Ofast -flto -fopenmp 
DEBUG_FLAGS = -g -O3 -DDEBUG 

# Default ENV_VARS
nt ?= 8
np ?= 1

# Create build
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Default target
all: serial parallel

# Release builds
build_serial: $(BUILD_DIR)/serial

$(BUILD_DIR)/serial: $(SRC_DIR)/serial.c $(INCLUDE_DIR)/serial.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(MATH) $(OPT_FLAGS)

build_parallel: $(BUILD_DIR)/parallel

$(BUILD_DIR)/parallel: $(SRC_DIR)/parallel.c $(INCLUDE_DIR)/parallel.h | $(BUILD_DIR)
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $< $(MATH) $(OPT_FLAGS)

# Debug builds
debug: serial_debug parallel_debug

serial_debug: $(BUILD_DIR)/serial_debug

$(BUILD_DIR)/serial_debug: $(SRC_DIR)/serial.c $(INCLUDE)/serial.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDES) -o $@ $< $(MATH)

parallel_debug: $(BUILD_DIR)/parallel_debug

$(BUILD_DIR)/parallel_debug: $(SRC_DIR)/parallel.c $(INCLUDE)/parallel.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDES) -o $@ $< $(MATH)

# Run
serial: build_serial
	$(BUILD_DIR)/serial -x 1000 -y 1000 -n 100

parallel: build_parallel
	export OMP_NUM_THREADS=$(nt); mpirun -np $(np) $(BUILD_DIR)/parallel -x 10000 -y 10000 -n 99

# Clean
clean: 
	rm -rf $(BUILD_DIR)

clean_all: clean
	rm -f *.o *.out core.*

