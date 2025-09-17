# =================================================================== 
# Makefile for compiling fftw_test.f90
# Adapted for Luke's cluster environment.
# ===================================================================

# Compiler: Using the MPI wrapper as in the main project
FC = mpif90

# Compiler flags from your project
FCFLAGS = -O2

# Include path for FFTW3 (taken from your Makefile)
INC = -I/software/fftw3/3.3.10/b1/include

# Library path and library for FFTW3 (taken from your Makefile)
LIB = -L/software/fftw3/3.3.10/b1/lib -lfftw3

# The target executable name
TARGET = fftw_test

# The source file
SOURCE = fftw_test.f90

# Default rule
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(SOURCE)
	$(FC) $(FCFLAGS) -o $(TARGET) $(SOURCE) $(INC) $(LIB)

# Rule to clean up the directory
.PHONY: clean
clean:
	rm -f $(TARGET)
