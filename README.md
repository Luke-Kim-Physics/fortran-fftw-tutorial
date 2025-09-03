# Fortran FFTW Tutorial

This repository contains a simple Fortran program that serves as a tutorial for using the FFTW3 library. The code demonstrates the Discrete Fourier Transform (DFT) of a 3D plane wave and numerically verifies the analytical result.

It is intended for students and researchers in computational physics, materials science, or any field that requires high-performance Fourier transforms in a compiled language.

---
## 📖 Background & Theory

This section provides a more detailed mathematical background on the concepts demonstrated in the code.

### 1. Plane Wave Representation: Continuous vs. Discrete

In solid-state physics, a plane wave is described by the continuous function $e^{i\vec{k} \cdot \vec{r}}$. In a computational setting, this function is discretized by sampling it on a finite grid within a simulation cell.

Let's consider a 3D simulation cell with dimensions $L_x, L_y, L_z$, which is discretized into an $N_x \times N_y \times N_z$ grid.

The **continuous position vector** is $\vec{r} = (x, y, z)$. Its discrete counterpart on the grid, indexed by $(i,j,k)$ where $i \in [0, N_x-1]$, etc., is:
```math
\vec{r}_{ijk} = \left( i\frac{L_x}{N_x}, \ j\frac{L_y}{N_y}, \ k\frac{L_z}{N_z} \right)
```

Due to periodic boundary conditions, the allowed **wavevectors ($\vec{k}$-vectors)** are quantized. A wavevector, indexed by integers $(q_x, q_y, q_z)$, is given by:
```math
\vec{k}_q = \left( q_x \frac{2\pi}{L_x}, \ q_y \frac{2\pi}{L_y}, \ q_z \frac{2\pi}{L_z} \right)
```
The phase of the plane wave, $\vec{k} \cdot \vec{r}$, for the discrete vectors is therefore:
```math
\vec{k}_q \cdot \vec{r}_{ijk} = 2\pi \left( \frac{q_x i}{N_x} + \frac{q_y j}{N_y} + \frac{q_z k}{N_z} \right)
```
This shows that the code's function, $f(i,j,k)$, is the exact discrete representation of the physical plane wave $\exp(i\vec{k}\_q \cdot \vec{r}\_{ijk})$.

### 2. The Discrete Fourier Transform of a Plane Wave

The forward 3D DFT, $\tilde{f}$, of a function $f(i,j,k)$ is defined as:
```math
\tilde{f}(n_x, n_y, n_z) = \sum_{i=0}^{N_x-1} \sum_{j=0}^{N_y-1} \sum_{k=0}^{N_z-1} f(i,j,k) e^{-2\pi i \left( \frac{n_x i}{N_x} + \frac{n_y j}{N_y} + \frac{n_z k}{N_z} \right)}
```
Substituting our discrete plane wave for $f(i,j,k)$, the expression becomes a product of three independent geometric series:
```math
\tilde{f}(n_x, n_y, n_z) = \left( \sum_{i=0}^{N_x-1} e^{2\pi i \frac{(q_x - n_x) i}{N_x}} \right) \times \left( \sum_{j=0}^{N_y-1} e^{2\pi i \frac{(q_y - n_y) j}{N_y}} \right) \times \left( \sum_{k=0}^{N_z-1} e^{2\pi i \frac{(q_z - n_z) k}{N_z}} \right)
```
Each of these sums follows the property:
```math
\sum_{i=0}^{N-1} e^{2\pi i \frac{m \cdot i}{N}} = 
\begin{cases} 
N & \text{if } m \text{ is a multiple of } N \text{ (i.e., } m=0 \text{ in this context)} \\
0 & \text{otherwise} 
\end{cases}
```
This can be expressed concisely using the Kronecker Delta, $\delta_{ab}$:
```math
\sum_{i=0}^{N_x-1} e^{2\pi i \frac{(q_x - n_x) i}{N_x}} = N_x \cdot \delta_{q_x, n_x}
```
Therefore, the final result for the DFT is a single non-zero point in reciprocal space:
```math
\tilde{f}(n_x, n_y, n_z) = (N_x \delta_{q_x, n_x}) (N_y \delta_{q_y, n_y}) (N_z \delta_{q_z, n_z}) = (N_x N_y N_z) \cdot \delta_{q_x, n_x} \delta_{q_y, n_y} \delta_{q_z, n_z}
```
This proves that the DFT of a single discrete plane wave is a single peak with a magnitude equal to the total number of grid points.

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
SOURCE = fftw_test.f90

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
