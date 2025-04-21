import java.util.Arrays;

public class Main {
    public static void main(String[] args) {
//        double[][] A = {
//                {1, -1, -4, 9},
//                {2, -3, 1, 5},
//                {1, 2, 0, -4},
//                {3, -2, -5, 1}
//        };
//        double[] B = {22, -3, -3, 3};

        double[][] A = {
                {4, -1, 0, 0},
                {-1, 3, -1, 0},
                {0, -1, 4, -1},
                {0, 0, -1, 3}
        };
        double[] B = {15, 10, 10, 10};

        System.out.println("Cramer:");
        double[] cramerSolution = solveByCramer(A, B);
        System.out.println(Arrays.toString(cramerSolution));

        System.out.println("Gauss-Seidel:");
        double[] gaussSeidelSolution = solveByGaussSeidel(A, B, 1e-6, 100);
        System.out.println(Arrays.toString(gaussSeidelSolution));
    }

    public static double[] solveByCramer(double[][] A, double[] B) {
        int n = A.length;
        double detA = determinant(A);
        if (detA == 0) throw new IllegalArgumentException("Determinant is zero, system has no unique solution.");

        double[] X = new double[n];
        for (int i = 0; i < n; i++) {
            double[][] Ai = new double[n][n];
            for (int r = 0; r < n; r++) {
                System.arraycopy(A[r], 0, Ai[r], 0, n);
                Ai[r][i] = B[r];
            }
            X[i] = determinant(Ai) / detA;
        }
        return X;
    }

    private static double determinant(double[][] matrix) {
        int n = matrix.length;
        if (n == 1) return matrix[0][0];
        if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

        double det = 0;
        for (int i = 0; i < n; i++) {
            double[][] minor = new double[n - 1][n - 1];
            for (int r = 1; r < n; r++) {
                for (int c = 0, colIndex = 0; c < n; c++) {
                    if (c == i) continue;
                    minor[r - 1][colIndex++] = matrix[r][c];
                }
            }
            det += matrix[0][i] * Math.pow(-1, i) * determinant(minor);
        }
        return det;
    }

    public static double[] solveByGaussSeidel(double[][] A, double[] B, double tol, int maxIter) {
        int n = A.length;
        double[] X = new double[n];
        double[] X_prev = new double[n];

        int k = 0;
        double maxDiff;

        do {
            maxDiff = 0;

            for (int i = 0; i < n; i++) {
                double sum = B[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= A[i][j] * X[j];
                    }
                }
                X_prev[i] = X[i];
                X[i] = sum / A[i][i];

                maxDiff = Math.max(maxDiff, Math.abs(X[i] - X_prev[i]));
            }

            k++;
        } while (maxDiff >= tol && k < maxIter);

        return X;
    }

}
