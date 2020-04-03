# Sparse Matrix Solver
Solves a randomly generated sparse matrix (that has an extremely high probability of being diagonally dominant) for a randomly generated array b, with size and density set by the user. 
Uses a stabilized Bi-conjugate gradient method for solving the system. This is an iterative approach. with the maximum number of iterations capped arbitrarily at 20. The solver runs twice, first without using a preconditioner, and a second time using a preconditioner equal to the diagonal of the generated matrix.

# Usage
Code is written in C#, and was compiled and run using Visual Studio.
