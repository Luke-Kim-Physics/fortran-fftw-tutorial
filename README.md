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
12312313213213213213213213213
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
This shows that the code's function, $f(i,j,k)$, is the exact discrete representation of the physical plane wave $\exp(i\vec{k}_q \cdot \vec{r}_{ijk})$.

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

\section{From Continuous Plane Waves to Discrete Fourier Transforms in Computational Physics}

\subsection{Plane Wave Representation: Continuous vs. Discrete}
In solid-state physics, a plane wave is described by the continuous function $e^{i\vec{k} \cdot \vec{r}}$. In a computational setting, this function is discretized by sampling it on a finite grid within a simulation cell.

\subsubsection{Vector Definitions}
Let's consider a 3D simulation cell with dimensions $L_x, L_y, L_z$, which is discretized into an $N_x \times N_y \times N_z$ grid.

The \textbf{continuous position vector} is $\vec{r} = (x, y, z)$. Its discrete counterpart on the grid, indexed by $(i,j,k)$ where $i \in [0, N_x-1]$, etc., is:
\begin{equation}
    \vec{r}_{ijk} = \left( i\frac{L_x}{N_x}, \ j\frac{L_y}{N_y}, \ k\frac{L_z}{N_z} \right)
\end{equation}

Due to periodic boundary conditions, the allowed \textbf{wavevectors (k-vectors)} are quantized. A wavevector, indexed by integers $(q_x, q_y, q_z)$, is given by:
\begin{equation}
    \vec{k}_q = \left( q_x \frac{2\pi}{L_x}, \ q_y \frac{2\pi}{L_y}, \ q_z \frac{2\pi}{L_z} \right)
\end{equation}

\subsubsection{The Dot Product}
The phase of the plane wave, $\vec{k} \cdot \vec{r}$, can now be calculated for the discrete vectors:
\begin{align*}
    \vec{k}_q \cdot \vec{r}_{ijk} &= \left( q_x \frac{2\pi}{L_x} \right) \left( i\frac{L_x}{N_x} \right) + \left( q_y \frac{2\pi}{L_y} \right) \left( j\frac{L_y}{N_y} \right) + \left( q_z \frac{2\pi}{L_z} \right) \left( k\frac{L_z}{N_z} \right) \\
    &= 2\pi \left( \frac{q_x i}{N_x} + \frac{q_y j}{N_y} + \frac{q_z k}{N_z} \right)
\end{align*}
This demonstrates that the discrete function used in the test code, $f(i,j,k) = \exp\left[ 2\pi i \left( \frac{q_x i}{N_x} + \dots \right) \right]$, is the exact discrete representation of the physical plane wave $\exp(i\vec{k}_q \cdot \vec{r}_{ijk})$.

\subsection{The Discrete Fourier Transform (DFT) of a Plane Wave}
The forward 3D DFT, $\tilde{f}$, of a function $f(i,j,k)$ is defined as:
\begin{equation}
    \tilde{f}(n_x, n_y, n_z) = \sum_{i=0}^{N_x-1} \sum_{j=0}^{N_y-1} \sum_{k=0}^{N_z-1} f(i,j,k) e^{-2\pi i \left( \frac{n_x i}{N_x} + \frac{n_y j}{N_y} + \frac{n_z k}{N_z} \right)}
\end{equation}
Substituting our discrete plane wave for $f(i,j,k)$:
\begin{equation}
    \tilde{f}(n_x, n_y, n_z) = \sum_{i,j,k} \exp\left[ 2\pi i \left( \frac{(q_x - n_x) i}{N_x} + \frac{(q_y - n_y) j}{N_y} + \frac{(q_z - n_z) k}{N_z} \right) \right]
\end{equation}
Using the property $e^{A+B+C}=e^A e^B e^C$, this triple summation separates into a product of three independent single sums:
\begin{equation}
    \tilde{f}(n_x, n_y, n_z) = S_x \cdot S_y \cdot S_z
\end{equation}
where, for the x-dimension:
\begin{equation}
    S_x = \sum_{i=0}^{N_x-1} \exp\left[ 2\pi i \frac{(q_x - n_x) i}{N_x} \right]
\end{equation}
This sum is a \textbf{geometric series} of the form $\sum_{i=0}^{N-1} r^i$ with a common ratio $r = \exp\left[ 2\pi i \frac{q_x - n_x}{N_x} \right]$.

\subsubsection{Analysis of the Geometric Series}
\begin{itemize}
    \item \textbf{Case 1: $q_x = n_x$}. The common ratio $r = \exp(0) = 1$. The sum becomes:
    \[ S_x = \sum_{i=0}^{N_x-1} (1)^i = N_x \]
    
    \item \textbf{Case 2: $q_x \neq n_x$}. The common ratio $r \neq 1$. Using the formula for the sum of a geometric series, $a(1-r^n)/(1-r)$:
    \[ S_x = \frac{1 - r^{N_x}}{1 - r} = \frac{1 - \left( \exp\left[ 2\pi i \frac{q_x - n_x}{N_x} \right] \right)^{N_x}}{1 - r} \]
    The numerator simplifies to:
    \[ 1 - \exp[2\pi i (q_x - n_x)] = 1 - 1 = 0 \]
    Since $e^{i*2\pi m} = cos(2\pi m) + isin(2\pi m) = 1+0 = 1$ (m = $q_x - n_x$, integer) and the denominator is non-zero, $S_x = 0$.
\end{itemize}
Combining these two cases using the Kronecker delta, $\delta_{ab}$:
\begin{equation}
    S_x = N_x \cdot \delta_{q_x, n_x}
\end{equation}
Applying this to all three dimensions, the final result for the DFT is:
\begin{equation}
    \tilde{f}(n_x, n_y, n_z) = (N_x \delta_{q_x, n_x}) (N_y \delta_{q_y, n_y}) (N_z \delta_{q_z, n_z}) = (N_x N_y N_z) \cdot \delta_{q_x, n_x} \delta_{q_y, n_y} \delta_{q_z, n_z}
\end{equation}
This proves that the DFT of a single discrete plane wave is a single non-zero point (a delta function) in reciprocal space, with a magnitude equal to the total number of grid points.

\subsection{Application to Bloch's Theorem}
The cell-periodic part of a Bloch function, $u_{\vec{k}}(\vec{r})$, can be expanded as a Fourier series over the reciprocal lattice vectors $\vec{G}$:
\begin{equation}
    u_{\vec{k}}(\vec{r}) = \sum_{\vec{G}} C_{\vec{k},\vec{G}} e^{i\vec{G}\cdot\vec{r}}
\end{equation}
The continuous Fourier coefficients, $C_{\vec{k},\vec{G}}$, are found by the inverse transform:
\begin{equation}
    C_{\vec{k},\vec{G}} = \frac{1}{\Omega} \int_{\text{cell}} u_{\vec{k}}(\vec{r}) e^{-i\vec{G}\cdot\vec{r}} d^3r
\end{equation}
where $\Omega$ is the cell volume.

To implement this numerically, the integral is replaced by a discrete sum over $N$ grid points:
\begin{align*}
    C_{\vec{k},\vec{G}} &\approx \frac{1}{\Omega} \sum_{j=0}^{N-1} u_{\vec{k}}(\vec{r}_j) e^{-i\vec{G}\cdot\vec{r}_j} \Delta V \quad \text{where } \Delta V = \frac{\Omega}{N} \\
    &\approx \frac{1}{\Omega} \left( \sum_{j=0}^{N-1} u_{\vec{k}}(\vec{r}_j) e^{-i\vec{G}\cdot\vec{r}_j} \right) \frac{\Omega}{N} \\
    &\approx \frac{1}{N} \sum_{j=0}^{N-1} u_{\vec{k}}(\vec{r}_j) e^{-i\vec{G}\cdot\vec{r}_j}
\end{align*}
The summation term is precisely the definition of the unnormalized forward DFT of the function $u_{\vec{k}}(\vec{r})$ evaluated at the grid points.
\begin{equation}
    C_{\vec{k},\vec{G}} \approx \frac{1}{N} \times \text{DFT}[u_{\vec{k}}(\vec{r}_j)]
\end{equation}
This confirms that to obtain the physically meaningful plane-wave coefficients, one must perform a forward DFT on the real-space function and divide by the total number of grid points, $N = n_x n_y n_z$.

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
