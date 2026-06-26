# Makefile for LetzElPhC Unit Tests
include ../../src/make.inc

# Re-use the include directories and flags from main make.inc
INCLUDES := -I../../src $(FFTW_INC) $(NETCDF_INC)
LIBS := $(LD_FLAGS) $(FFTW3_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(HDF5_LIB) -lm

SRCS := $(wildcard *.c)
BINS := $(SRCS:.c=)

# We also need to link against the object files in the src directory 
# Assuming they are already built when running unit tests
SRC_OBJS := $(wildcard ../../src/*/*.o ../../src/*.o)
# Filter out main.o if it exists so we don't get multiple main() definitions
SRC_OBJS := $(filter-out ../../src/lelphc.o ../../src/main.o, $(SRC_OBJS))

all: $(BINS)

%: %.c
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) $(INCLUDES) $< $(SRC_OBJS) -o $@ $(LIBS)

clean:
	rm -f $(BINS)
