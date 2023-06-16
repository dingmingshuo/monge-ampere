# Monge-Ampére Equation Solver

## Introduction


In this repository, we are going to implement two-scale method for Monge-Ampére equation on multi-thread CPU.
The equation is:
```math
\left\{
\begin{aligned}
    \det D^2 u = f(x) &\qquad\text{in} &\Omega \subset \mathbb{R}^2, \\
    u = g(x) &\qquad\text{on} &\partial \Omega
\end{aligned}
\right.
```
where $\Omega$ is a domain in $\mathbb{R}^2$ and $f(x), g(x)$ are given functions.
In this repository, we assume that $\Omega$ is a square $[0, 1]^2$, and discretized by a unstructured grid with diagram equals to $h$.


## Compile and Run

### Install and Compile

This code uses [Eigen](https://eigen.tuxfamily.org/index.php) as the linear algebra library. To compile the code, you need to install Eigen **first**. To install Eigen, run the following command in the terminal:

```bash
git submodule update --init --recursive
```

To compile the code and run the tests, run the following command in the terminal:

```bash
make all
```

### Run Examples

If you want to run the examples without compiling all tests, run the following command in the terminal:

```bash
make ./build/example/perron_example_1 && ./build/example/perron_example_1
```

Following examples is available:

- ./build/example/perron_example_1
- ./build/example/perron_example_2
- ./build/example/perron_example_3
- ./build/example/newton_example_1
- ./build/example/newton_example_2
- ./build/example/newton_example_3

Super parameters can be changed in the corresponding source `.cpp` file.

## Code Structure

```
.
├── example
│   ├── newton_example_1.cpp    // Example 1 for Newton method
│   ├── newton_example_2.cpp    // Example 2 for Newton method
│   ├── newton_example_3.cpp    // Example 3 for Newton method
│   ├── perron_example_1.cpp    // Example 1 for Perron iteration
│   ├── perron_example_2.cpp    // Example 2 for Perron iteration
│   └── perron_example_3.cpp    // Example 3 for Perron iteration
├── include
│   └── ma.hpp                  // Header file
├── Makefile
├── mesh                        // Mesh files (Typo: 2e-k means h = 2^{-k})
│   ├── 2e-2.obj
│   ├── 2e-3.obj
│   ├── 2e-4.obj
│   ├── 2e-5.obj
│   ├── 2e-6.obj
│   ├── 2e-7.obj
│   └── 2e-8.obj
├── README.md
├── source
│   ├── mesh.cpp                // Mesh struct
│   ├── mesh_function.cpp       // MeshFunction struct
│   ├── point.cpp               // Point struct
│   ├── solver.cpp              // Solver namespace (include CG and GMRES)
│   └── two_scale.cpp           // TwoScale namespace (include Newton and Perron)
└── test
    ├── mesh.cpp                // Test for Mesh struct
    ├── mesh_function.cpp       // Test for MeshFunction struct
    ├── solver.cpp              // Test for Solver namespace
    └── test_main.cpp           // Main function for all tests
```
