using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace SparseMatrixSolver
{
    
    class PreconditionedMatrixSolver
    {
        //static void Main(string[] args)
        //{
        //    Console.WriteLine("Hello World!");
        //    int size = 5000;
        //    var matrix = GenerateMatrix(size, 1);
        //    var b = new double[size];
        //    for(int i = 0; i < size; i ++)
        //    {
        //        b[i] = i + 1;
        //    }
        //    var x = PreconditionedBiCGSTABSolver(size, b, matrix, 20, 0.01);
        //    //for (int i = 0; i < size; i++)
        //    //{
        //    //    Console.WriteLine(x[i]);
        //    //}
        //    Console.WriteLine(ValidateSolution(matrix, x, b));
        //}

        public static List<MatrixEntry> GenerateMatrix(int size, int density)
        {
            var matrix = new List<MatrixEntry>();
            var random = new Random();
            //for (int i = 0; i < size; i++)
            //{
            //    for(int j = 0; j < size; j++)
            //    {
            //        if(i != j)
            //        {
            //            var randomNumber = random.Next(1, 1000);
            //            if (randomNumber <= 10 * density)
            //            {
            //                matrix.Add(new MatrixEntry(randomNumber, i, j));
            //            }
            //        }
            //    }
            //}
            for (int i = 0; i < size; i++)
            {
                matrix.Add(new MatrixEntry(random.Next(200 * density, 300 * density), i, i));
            }
            for (int i = 1; i < size; i++)
            {
                matrix.Add(new MatrixEntry(random.Next(50 * density, 60 * density), i, i - 1));
            }
            for (int i = 0; i < size - 1; i++)
            {
                matrix.Add(new MatrixEntry(random.Next(50 * density, 60 * density), i, i + 1));
            }
            matrix.Add(new MatrixEntry(1, 0, 1));
            return matrix;
        }

        public static double[] PreconditionedBiCGSTABSolver(int size, double[] b, List<MatrixEntry> A, int maxIterations, double precision)
        {
            var xPrevious = new double[size];
            var r0 = ArraySubtraction(b, CompressedArrayMultiply(A, xPrevious));
            var r = new double[size];
            for (int i = 0; i < size; i++)
            {
                r[i] = r0[i];
            }
            var alpha = 1.0;
            var rhoPrevious = 1.0;
            var omegaPrevious = 1.0;
            var vPrevious = new double[size];
            var pPrevious = new double[size];

            for (int i = 0; i < maxIterations; i++)
            {
                var rhoCurrent = ArrayDotProduct(r0, r);
                var beta = (rhoCurrent / rhoPrevious) / (alpha / omegaPrevious);
                var dummy = ArraySubtraction(pPrevious, ArrayScalarMultiplication(omegaPrevious, vPrevious));
                var p = ArrayAddition(r, ArrayScalarMultiplication(beta, dummy));
                var v = CompressedArrayMultiply(A, p);
                alpha = rhoCurrent / (ArrayDotProduct(r0, v));
                var h = ArrayAddition(xPrevious, ArrayScalarMultiplication(alpha, p));
                if (ValidateSolution(A, h, b) < precision)
                {
                    return h;
                }

                var s = ArraySubtraction(r, ArrayScalarMultiplication(alpha, v));
                var t = CompressedArrayMultiply(A, s);
                var omega = ArrayDotProduct(t, s) / ArrayDotProduct(t, t);

                var x = ArrayAddition(h, ArrayScalarMultiplication(omega, s));
                Console.WriteLine("Error at iteration " + i + ": " + ValidateSolution(A, x, b));

                if (ValidateSolution(A, x, b) < precision)
                {
                    return x;
                }
                r = ArraySubtraction(s, ArrayScalarMultiplication(omega, t));
                xPrevious = x;
                rhoPrevious = rhoCurrent;
                omegaPrevious = omega;
                vPrevious = v;
                pPrevious = p;
            }

            return xPrevious;
        }

        public static double ArrayDotProduct(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for dot product");
            }
            return array1.Zip(array2, (d1, d2) => d1 * d2)
                        .Sum();
        }

        public static double[] ArraySubtraction(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for subtraction");
            }
            return array1.Zip(array2, (d1, d2) => d1 - d2).ToArray();
        }

        public static double[] ArrayAddition(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
            {
                throw new Exception("Arrays must be the same length for subtraction");
            }
            return array1.Zip(array2, (d1, d2) => d1 + d2).ToArray();
        }

        public static double[] ArrayScalarMultiplication(double scalar, double[] array)
        {
            return array.Select(value => scalar * value).ToArray();
        }

        public static double[] CompressedArrayMultiply(List<MatrixEntry> A, double[] x)
        {
            var result = new double[x.Length];
            foreach (var matrixEntry in A)
            {
                result[matrixEntry.RowIndex] += matrixEntry.Value * x[matrixEntry.ColumnIndex];
            }
            return result;
        }

        public static double ValidateSolution(List<MatrixEntry> A, double[] x, double[] b)
        {
            var result = CompressedArrayMultiply(A, x);
            double residualSquared  = 0;
            for(int i = 0; i < x.Length; i++)
            {
                residualSquared += Math.Pow((result[i] - b[i]), 2);
            }
            return residualSquared;
        }
    }
}
