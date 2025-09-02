# Fortran FFTW Tutorial

This repository contains a simple Fortran program that serves as a tutorial for using the FFTW3 library. The code demonstrates the Discrete Fourier Transform (DFT) of a 3D plane wave and numerically verifies the analytical result.

It is intended for students and researchers in computational physics, materials science, or any field that requires high-performance Fourier transforms in a compiled language.

---

## 📖 Background & Theory

In physics, a plane wave is described by the continuous function $e^{i\vec{k} \cdot \vec{r}}$. When we discretize this function on a finite grid, its Discrete Fourier Transform is expected to be a single, non-zero point (a Kronecker delta function) in reciprocal space.

This code verifies this fundamental principle. It:
1.  **Creates** a complex plane wave function on a 3D grid in real space.
2.  **Transforms** this function to reciprocal space using a forward 3D FFT call from the FFTW3 library.
3.  **Verifies** that the result is a single peak at the correct frequency, with a magnitude equal to the total number of grid points ($N = N_x \times N_y \times N_z$).
4.  **Transforms** the result back to real space using a backward FFT to recover the original function, confirming the accuracy of the round-trip process.

---

## ⚙️ Requirements

To compile and run this code, you will need:
* A Fortran compiler (e.g., `gfortran`, `ifort`, `mpif90`).
* The **FFTW3 library** installed on your system. You must know the path to its include and lib directories.

---

## 🚀 How to Compile and Run

### 1. The Makefile
This repository includes a `Makefile` for easy compilation. Before running `make`, you may need to edit the `INC` and `LIB` variables inside the `Makefile` to match the paths where FFTW3 is installed on your system.

```makefile
# Makefile for fftw_verification_test.f90
FC = gfortran
FCFLAGS = -O2

# --- EDIT THESE PATHS FOR YOUR SYSTEM ---
INC = -I/path/to/your/fftw3/include
LIB = -L/path/to/your/fftw3/lib -lfftw3
# ----------------------------------------

TARGET = fftw_test
SOURCE = fftw_verification_test.f90

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(FC) $(FCFLAGS) -o $(TARGET) $(SOURCE) $(INC) $(LIB)

clean:
	rm -f $(TARGET)
```

### 2. Compilation
Once the `Makefile` is configured, simply run the `make` command in your terminal:
```bash
make
```
This will create an executable file named `fftw_test`.

### 3. Execution
Run the compiled program:
```bash
./fftw_test
```

---

## 📊 Expected Output & Explanation

Running the program should produce the following output. Here is an explanation of what each part means.

```
 --- Verification of Forward Fourier Transform ---
 Grid size:            8x          8x          8 =         512 points.
 Input plane wave frequency indices (qx,qy,qz):            2           3           4
 Theoretical peak value should be close to:    512.00000000000000

 Peak found at index (  3  4  5) with value:   512.000000     -0.000000
 Value at index      (  1  1  1) is:                 -0.000000     -0.000000
 Value at index      (  2  2  2) is:                 -0.000000      0.000000

 --- Verification of Backward Fourier Transform ---
 Comparing original function with the back-transformed one.
 The difference should be close to zero.
 Original value at (5,5,5):           -1.000000      0.000000
 Back-transformed value at (5,5,5):   -1.000000      0.000000
 Absolute difference:                 0.1972E-30
```

* **Peak Position**: The input frequency indices were `(2,3,4)`. Since Fortran arrays are 1-based, the peak correctly appears at index `(3,4,5)`.
* **Peak Value**: The magnitude of the peak is `512`, which is exactly the total number of grid points ($8 \times 8 \times 8$), as predicted by theory. Other points are numerically zero.
* **Backward Transform**: The program successfully recovers the original function. The "Absolute difference" of $\approx 10^{-31}$ is effectively zero and represents the limit of machine precision, confirming a perfect round-trip transformation.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
## Author
Luke Kim
