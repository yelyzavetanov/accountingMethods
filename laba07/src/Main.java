public class Main {

    // Функція для знаходження значення факторіала
    static int factorial(int n) {
        int fact = 1;
        for (int i = 2; i <= n; i++) {
            fact *= i;
        }
        return fact;
    }

    // Функція для побудови таблиці різниць
    static void computeForwardDifference(double[][] differenceTable, double[] y, int n) {
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < n - i; j++) {
                differenceTable[j][i] = differenceTable[j + 1][i - 1] - differenceTable[j][i - 1];
            }
        }
    }

    // Функція для інтерполяції (формула Ньютона вперед)
    static double interpolate(double value, double[] x, double[][] differenceTable, int n) {
        double h = x[1] - x[0];  // крок
        double p = (value - x[0]) / h;
        double result = differenceTable[0][0];

        for (int i = 1; i < n; i++) {
            double term = differenceTable[0][i];
            for (int j = 0; j < i; j++) {
                term *= (p - j);
            }
            term /= factorial(i);
            result += term;
        }
        return result;
    }

    public static void main(String[] args) {
        // Відомі вузли
        double[] x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double[] y = {-1, 2, 17, 50, 107, 194, 317, 482, 695, 962, 1289};
        int n = x.length;

        // Створюємо таблицю різниць
        double[][] differenceTable = new double[n][n];

        // Ініціалізуємо перший стовпець таблиці значеннями y
        for (int i = 0; i < n; i++) {
            differenceTable[i][0] = y[i];
        }

        // Обчислюємо таблицю різниць
        computeForwardDifference(differenceTable, y, n);

        // Точки, де потрібно обчислити значення
        double[] points = {1.5, 3.4, 2.8};

        System.out.println("Результати інтерполяції:");
        for (double point : points) {
            double result = interpolate(point, x, differenceTable, n);
            System.out.printf("f(%.1f) ≈ %.4f%n", point, result);
        }
    }
}
