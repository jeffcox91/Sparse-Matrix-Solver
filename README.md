# Sparse Matrix Solver
Solves a randomly generated sparse matrix (with all diagonal entries guaranteed to be full, and much larger than other entries in the matrix) for a randomly generated array b, with size and density set by the user. 
Uses a stabilized bi-conjugate gradient method for solving the system. This is an iterative approach. with the maximum number of iterations capped arbitrarily at 20. The solver runs twice, first without using a preconditioner, and a second time using a preconditioner equal to the diagonal of the generated matrix.

More information about this method can be found in the paper "Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems" by H. A. van der Vorst, which can be found at https://epubs.siam.org/doi/10.1137/0913035. Publicly available information about the method can also be found at https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method#endnote_bicgstab2

# Usage
Code is written in C#, and was compiled and run using Visual Studio.
