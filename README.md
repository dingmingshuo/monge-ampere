# Monge-Ampére Equation Solver

## Introduction


In this repository, we are going to implement two-scale method for Monge-Ampére equation on multi-thread CPU.
The equation is:
$$
\left\{
\begin{aligned}
    \det D^2 u = f(x) &\qquad\text{in} &\Omega \subset \mathbb{R}^2, \\
    u = g(x) &\qquad\text{on} &\partial \Omega
\end{aligned}
\right.
$$
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

Super parameters can be changed in the corresponding source `.cpp` file.

## Code Structure