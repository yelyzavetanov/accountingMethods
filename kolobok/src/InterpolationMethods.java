import java.util.Arrays;

public class InterpolationMethods {

    // Різницевий метод за першою формулою Ньютона (рівновіддалені вузли)
    public static double newtonFirst(double[] x, double[] y, double xVal) {
        int n = x.length;
        double[][] diff = new double[n][n];

        // Копіюємо перший стовпець
        for (int i = 0; i < n; i++) {
            diff[i][0] = y[i];
        }

        // Обчислення різниць
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
            }
        }

        double h = x[1] - x[0];
        double t = (xVal - x[0]) / h;
        double result = y[0];
        double term = 1;

        for (int i = 1; i < n; i++) {
            term *= (t - (i - 1)) / i;
            result += term * diff[0][i];
        }

        return result;
    }

    // Різницевий метод за другою формулою Ньютона (рівновіддалені вузли)
    public static double newtonSecond(double[] x, double[] y, double xVal) {
        int n = x.length;
        double[][] diff = new double[n][n];

        for (int i = 0; i < n; i++) {
            diff[i][0] = y[i];
        }

        for (int j = 1; j < n; j++) {
            for (int i = n - 1; i >= j; i--) {
                diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
            }
        }

        double h = x[1] - x[0];
        double t = (xVal - x[n - 1]) / h;
        double result = y[n - 1];
        double term = 1;

        for (int i = 1; i < n; i++) {
            term *= (t + (i - 1)) / i;
            result += term * diff[n - 1][i];
        }

        return result;
    }

    // Метод Лагранжа (працює з будь-якими вузлами)
    public static double lagrange(double[] x, double[] y, double xVal) {
        int n = x.length;
        double result = 0;

        for (int i = 0; i < n; i++) {
            double term = y[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    term *= (xVal - x[j]) / (x[i] - x[j]);
                }
            }
            result += term;
        }

        return result;
    }

    public static void main(String[] args) {
        double[] x = {1, 2, 3, 4};
        double[] y = {1, 4, 9, 16}; // y = x^2
        double xVal = 2.5;

        System.out.printf("Ньютон (1-а формула): %.6f\n", newtonFirst(x, y, xVal));
        System.out.printf("Ньютон (2-а формула): %.6f\n", newtonSecond(x, y, xVal));
        System.out.printf("Лагранж: %.6f\n", lagrange(x, y, xVal));
    }
}
