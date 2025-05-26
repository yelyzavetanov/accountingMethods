import java.io.*;
import java.util.*;

public class Main {

    static final int N = 4;
    static final double EPS = 0.001;

    public static void main(String[] args) {
        double[][] A = new double[N][N];
        double[] b = new double[N];

        System.out.println("Поточна директорія: " + new File(".").getAbsolutePath() + "\n");

        try {
            readMatrixFromFile( "C:\\comp\\vntu\\accountingMethods\\kursach\\src\\matrix.txt", A, b);
        } catch (IOException e) {
            System.err.println("Помилка при читанні файлу: " + e.getMessage());
            return;
        }

        double[][] B = {
                {0.2, -0.1},
                {0.1,  0.3}
        };
        System.out.println("Норма 1: " + norm1(A));
        System.out.println("Норма 2: " + norm2(A));
        System.out.println("E-Норма: " + eNorm(A));


        System.out.println("Гауссa-Жорданa:");
        double[] gaussJordanSolution = gaussJordan(A, b);
        printSolution(gaussJordanSolution);

        System.out.println("\nІтераційний метод:");
        double[] iterativeSolution = iterativeMethod(A, b);
        printSolution(iterativeSolution);
    }

    public static void readMatrixFromFile(String filename, double[][] A, double[] b) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line;
        int row = 0;

        while ((line = reader.readLine()) != null && row < N) {
            String[] parts = line.trim().split("\\s+");
            for (int i = 0; i < N; i++) {
                A[row][i] = Double.parseDouble(parts[i]);
            }
            b[row] = Double.parseDouble(parts[N]);
            row++;
        }

        reader.close();
    }

    public static double[] gaussJordan(double[][] A, double[] b) {
        int n = A.length;
        double[][] M = new double[n][n + 1];

        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = b[i];
        }

        for (int i = 0; i < n; i++) {
            double div = M[i][i];
            for (int j = 0; j <= n; j++) {
                M[i][j] /= div;
            }

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = M[k][i];
                    for (int j = 0; j <= n; j++) {
                        M[k][j] -= factor * M[i][j];
                    }
                }
            }
        }

        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = M[i][n];
        }
        return x;
    }

    public static double[] iterativeMethod(double[][] A, double[] b) {
        double[] x = new double[N];
        double[] prevX = new double[N];
        boolean converged;
        int iterations = 0;

        do {
            iterations++;
            System.arraycopy(x, 0, prevX, 0, N);
            for (int i = 0; i < N; i++) {
                double sum = b[i];
                for (int j = 0; j < N; j++) {
                    if (i != j) {
                        sum -= A[i][j] * prevX[j];
                    }
                }
                x[i] = sum / A[i][i];
            }

            converged = true;
            for (int i = 0; i < N; i++) {
                if (Math.abs(x[i] - prevX[i]) > EPS) {
                    converged = false;
                    break;
                }
            }

        } while (!converged);

        System.out.println("Кількість ітерацій: " + iterations);
        return x;
    }

    public static double norm1(double[][] B) {
        double maxSum = 0;
        for (double[] row : B) {
            double rowSum = 0;
            for (double val : row) {
                rowSum += Math.abs(val);
            }
            maxSum = Math.max(maxSum, rowSum);
        }
        return maxSum;
    }

    public static double norm2(double[][] B) {
        int cols = B[0].length;
        double maxSum = 0;
        for (int j = 0; j < cols; j++) {
            double colSum = 0;
            for (double[] row : B) {
                colSum += Math.abs(row[j]);
            }
            maxSum = Math.max(maxSum, colSum);
        }
        return maxSum;
    }

    public static double eNorm(double[][] B) {
        double sum = 0;
        for (double[] row : B) {
            for (double val : row) {
                sum += val * val;
            }
        }
        return Math.sqrt(sum);
    }

    public static void printSolution(double[] x) {
        for (int i = 0; i < x.length; i++) {
            System.out.printf("x%d = %.5f\n", i + 1, x[i]);
        }
    }
}