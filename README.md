<!-- ABOUT THE PROJECT -->
## About The Project

This is an implementation of Fast Multipole Method in C that supports `multi-precision` floating point computation for solving the following question.
<br>
<br>
Given Sources points $s_i$ for $i = 1, ..., m$, Target points $t_j$ for $j = 1,..., n$, and scalar $a_i$ for $i = 1, ..., m$ compute $F(t_j) = \sum_{i = 1}^m a_i\cdot K(t_j, s_i)$ for $j = 1,..., m$ where:

* $K(y, x) = \frac{1}{y - x}$ or
* $K(y, x) = \log |y - x|$.

Naive computation has complexity $O(mn)$ and Fast Multipole Method can compute approximations in $O(m) + O(n)$. This implementation adopts an adaptive partitioning approach to accommodate all distributions for the source points.


## Built With

* [A Short Course on Fast Multipole Methods](https://math.nyu.edu/~greengar/shortcourse_fmm.pdf)
* [GMP](https://gmplib.org/)
* [MPFR](https://www.mpfr.org/)
* [MPC](https://www.multiprecision.org/mpc/)

## How To Install

The users would need to compile the program themselves, and
this program relies on packages GMP, MPFR, and MPC and will not compile without them. 
For example with GCC compiler 4.2.1 the program can be compile from the root directory with  
<br>
`gcc main.c src/aux.c src/const.c src/expansion.c src/treeMP.c -lm -lgmp -lmpfr -lmpc`
<br>

## How To Use

The program reads inputs $s_i, t_i$, and $a_i$ from files, does the computation, and writes the output to file.  
The executable take the following required arguments

*    `-s filename`:        set filename for the input of sources
*    `-t filename`:        set filename for the input of targets
*    `-o filename`:        set filename to store the output. 

and following optional arguments

*    `-v filename`:        set filename for the input of scalars. 
                       All scaler are set to $1$ if this file is not supplied.
*    `-m binary`:          set the kernel function with 
                       0 for $ K(y, x) = 1/(y-x)$ and 1 for $K(y, x) = \log (y - x)$.
*    `-d integer`:          set the degree (uint) of multipole expansion. 
                       The default and minimum is 15 (and maximum is 25 for now).
*    `-l integer`:          set the degree (uint) of local expansion. 
                       The default and minimum is 20 (and maximum is 30 for now).
*    `-x integer`:          set the maximum source number in the finest grids, 
                       i.e. a grid will not be further sub-divided if number of sources 
                       contained is less than this number. The default and minimum is 300.
*    `-p integer`:          set the number of digits used in computation. Double precision used 53 digits.  

For example
<br>
`/a.out -s examples/easy10000 -t examples/target5000 -o res_easy_dp.txt -d 20 -l 25 -x 50 -p 256`
<br>


## Input File Format

For computation with double precision $(p \le 53)$, the double complex numbers are read
with `scanf(..., "(%lf %lf)",...)`.
<br>
For computation with multi-precision, the multi-precision complex numbers are read with `mpc_inp_str` from MPC. 
<br>
Please check `/examples/easy10000` and `/examples/uniform10000` for samples input files.

## Future Plan

* Extend the degree bound of Multipole Expansion and Local Expansion.
* Implement parallelizing using `pthread`.

