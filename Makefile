# Compiler
CC = gcc
MPICC = mpicc
ARCH ?= native

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# Include paths
INCLUDES = -I$(INCLUDE_DIR)

# Flags
MATH = -lm
CFLAGS = -std=c99 -Wall -Wextra -Wno-unknown-pragmas
OPT_FLAGS = -O3 -fopenmp -march=${ARCH}

# Default ENV_VARS
nt ?= 8
np ?= 1
x ?= 10000
y ?= 10000
p ?= 0
n ?= 100
f ?= 10
e ?= 4 
E ?= 1.0
o ?= 0
v ?= 0

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

# Run
serial: build_serial
	$(BUILD_DIR)/serial -x $(x) -y $(y) -n $(n) -p $(p) -e $(e) -E $(E) -f $(f) -o $(o) 

parallel: build_parallel
	export OMP_NUM_THREADS=$(nt); mpirun -np $(np) $(BUILD_DIR)/parallel -x $(x) -y $(y) -n $(n) -p $(p) -e $(e) -E $(E) -f $(f) -o $(o) -v $(v)

# Clean
clean: 
	rm -rf $(BUILD_DIR)

clean_all: clean
	rm -f *.o *.out core.*

