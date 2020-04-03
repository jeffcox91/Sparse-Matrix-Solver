using System;
using static SparseMatrixSolver.MatrixFunctions;

namespace SparseMatrixSolver
{
    class SparseMatrixSolver
    {
        static void Main(string[] args)
        {
            var size = 0;
            var density = 0.0;
            while (true)
            {
                Console.WriteLine("Enter size of matrix to solve (up to 100000, -1 to exit program):");
                size = Convert.ToInt32(Console.ReadLine());
                if (size < 0)
                {
                    break;
                }

                Console.WriteLine("Enter approximate desired density of matrix (0-100):");
                density = Convert.ToDouble(Console.ReadLine());
                if(density < 0 || density > 100)
                {

                }
                
                var diagonal = GenerateDiagonalEntries(size);
                var matrix = GenerateMatrix(size, density, diagonal);

                var b = new double[size];
                var random = new Random();
                for (int i = 0; i < size; i++)
                {
                    b[i] = random.Next(1, 100);
                }

                var x = BiCGSTABSolver(size, b, matrix, 20, Math.Pow(10, -11));
                Console.WriteLine("Final error with no preconditioner: " + SumSquaredErrors(matrix, x, b));
                Console.WriteLine("Alternate approach with preconditioner: ");
                var x_alt = BiCGSTABSolver(size, b, matrix, 20, Math.Pow(10, -11), diagonal);
                Console.WriteLine("Final error with preconditioner: " + SumSquaredErrors(matrix, x_alt, b));
            }
        }
    }
}
