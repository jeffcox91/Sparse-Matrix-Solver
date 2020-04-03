using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SparseMatrixSolver
{
    struct MatrixEntry
    {
        public double Value;
        public int RowIndex;
        public int ColumnIndex;

        public MatrixEntry(double value, int rowIndex, int columnIndex)
        {
            Value = value;
            RowIndex = rowIndex;
            ColumnIndex = columnIndex;
        }

    }
    static class MatrixFunctions
    {
        /*
         * Produces a random matrix, that has an extremely high probability of being diagonally dominant.
         * The diagonal entries are included as an optional argument, in case they have been pre-calculated to be used as a pre-conditioner.
         */
        public static List<MatrixEntry> GenerateMatrix(int size, double density, double[] diagonalEntries = null)
        {
            var matrix = new List<MatrixEntry>();
            var random = new Random();
            var densityAsPercentage = density / 100;

            /*
             * Loops over all possible entries in the matrix and generates a random number between 0 and 1.
             * If that is less than the density of the matrix, adds a new Matrix Entry with a value between 0 and 10.
             */
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (i != j)
                    {
                        var randomNumber = random.NextDouble();
                        if (randomNumber <= densityAsPercentage)
                        {
                            matrix.Add(new MatrixEntry(randomNumber * 10, i, j));
                        }
                    }
                }
            }

            if (diagonalEntries == null)
            {
                diagonalEntries = GenerateDiagonalEntries(size);
            }

            for (int i = 0; i < size; i++)
            {
                matrix.Add(new MatrixEntry(diagonalEntries[i], i, i));
            }

            return matrix;
        }

        /*
         * This is separated out to obtain the preconditioner, which in the case of a diagonally dominant matrix is simply the diagonal of the matrix
         */
        public static double[] GenerateDiagonalEntries(int size)
        {
            var diagonalEntries = new double[size];
            var random = new Random();
            for (int i = 0; i < size; i++)
            {
                diagonalEntries[i] = (random.NextDouble() + 1) * size; //For large, sparse matrices, virtually guarantees a diagonally dominant matrix
            }
            return diagonalEntries;
        }

        /*
         * Implements a stabilized bi-conjugate gradient solver, as described in more detail here: https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method#endnote_bicgstab2
         * 
         */
        public static double[] BiCGSTABSolver(int size, double[] b, List<MatrixEntry> A, int maxIterations, double precision, double[] preconditioner = null)
        {
            var xPrevious = new double[size]; //Initial guess.
            var r = ArraySubtraction(b, CompressedArrayMultiply(A, xPrevious)); //Initial guess for the residual. 
            var r_hat = new double[size];
            var p = new double[size];
            var v = new double[size];
            for (int i = 0; i < size; i++)
            {
                r_hat[i] = r[i]; // r_hat is an arbitrary vector, but the dot product of r and r_hat must be non-zero, so a simple choice for r_hat is to set it equal to r0.
            }
            var alpha = 1.0;
            var rhoPrevious = 1.0;
            var omega = 1.0;

            for (int i = 0; i < maxIterations; i++)
            {
                var rhoCurrent = ArrayDotProduct(r_hat, r); //Direction vector
                var beta = (rhoCurrent / rhoPrevious) / (alpha / omega);

                var temp = ArraySubtraction(p, ArrayScalarMultiplication(omega, v));
                p = ArrayAddition(r, ArrayScalarMultiplication(beta, temp));

                double[] conditioned_p;
                if(preconditioner == null)
                {
                    conditioned_p = p;
                }
                else
                {
                    conditioned_p = ArrayElementDivide(p, preconditioner);
                }
                v = CompressedArrayMultiply(A, conditioned_p);
                alpha = rhoCurrent / (ArrayDotProduct(r_hat, v));

                //Early convergence condition.
                var h = ArrayAddition(xPrevious, ArrayScalarMultiplication(alpha, conditioned_p));
                if (SumSquaredErrors(A, h, b) < precision)
                {
                    return h;
                }

                //s & t act as stabilizers
                var s = ArraySubtraction(r, ArrayScalarMultiplication(alpha, v));

                double[] conditioned_s;
                if (preconditioner == null)
                {
                    conditioned_s = s;
                }
                else
                {
                    conditioned_s = ArrayElementDivide(s, preconditioner);
                }
                var t = CompressedArrayMultiply(A, conditioned_s);
                omega = ArrayDotProduct(t, s) / ArrayDotProduct(t, t);

                var x = ArrayAddition(h, ArrayScalarMultiplication(omega, conditioned_s));
                Console.WriteLine("Error at iteration " + i + ": " + SumSquaredErrors(A, x, b));

                if (SumSquaredErrors(A, x, b) < precision)
                {
                    return x;
                }
                r = ArraySubtraction(s, ArrayScalarMultiplication(omega, t));
                xPrevious = x;
                rhoPrevious = rhoCurrent;
            }

            return xPrevious;
        }


        public static double ArrayDotProduct(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for dot product");
            }
            return array1.Zip(array2, (element1, element2) => element1 * element2)
                        .Sum();
        }

        /*
         * Subtracts each element from array 2 from the corresponding element in array 1. Arrays must be the same length
         */
        public static double[] ArraySubtraction(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for subtraction");
            }
            return array1.Zip(array2, (element1, element2) => element1 - element2).ToArray();
        }

        /*
         * Adds two arrays together, element by element. Arrays must be the same length
         */
        public static double[] ArrayAddition(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for addition");
            }
            return array1.Zip(array2, (element1, element2) => element1 + element2).ToArray();
        }

        /*
         * Multiplies each element of an array by a scalar
         */
        public static double[] ArrayScalarMultiplication(double scalar, double[] array)
        {
            return array.Select(value => scalar * value).ToArray();
        }

        /*
         * Divides each element of array 1 by the corresponding element in array 2. Arrays must be the same length
         */
        public static double[] ArrayElementDivide(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for element-wise divide");
            }
            return array1.Zip(array2, (element1, element2) => element1 / element2).ToArray();
        }

        /*
         * This function performs a standard array multiplication operation, but using the compressed array format.
         * Due to the limitations of the compressed array format, this function has no built-in checks to verify that there are no column indexes
         * in the compressed array that are greater than the size of the array x.
         */
        public static double[] CompressedArrayMultiply(List<MatrixEntry> A, double[] x)
        {
            var result = new double[x.Length];
            foreach (var matrixEntry in A)
            {
                result[matrixEntry.RowIndex] += matrixEntry.Value * x[matrixEntry.ColumnIndex];
            }
            return result;
        }

        /*
         * Calculates the sum of the squared differences between the result of A*x and the actual answer b
         */
        public static double SumSquaredErrors(List<MatrixEntry> A, double[] x, double[] b)
        {
            var result = CompressedArrayMultiply(A, x);
            double residualSquared = 0;
            for (int i = 0; i < x.Length; i++)
            {
                residualSquared += Math.Pow((result[i] - b[i]), 2);
            }
            return residualSquared;
        }
    }
}
