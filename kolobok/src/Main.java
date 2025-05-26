public class Main {

    public static void main(String[] args) {
        // Тестування методів
        double[][] A = {
                {1, 2},
                {3, 4}
        };
        double[][] B = {
                {5, 6},
                {7, 8}
        };

        System.out.println("а) Складання матриць:");
        printMatrix(addMatrices(A, B));

        System.out.println("а) Віднімання матриць:");
        printMatrix(subtractMatrices(A, B));

        System.out.println("б) Множення матриці на коефіцієнт:");
        printMatrix(multiplyMatrixByScalar(A, 2));

        System.out.println("в) Множення матриць:");
        printMatrix(multiplyMatrices(A, B));

        System.out.println("г) Транспонування матриці:");
        printMatrix(transposeMatrix(A));

        System.out.println("д) Обернена матриця:");
        printMatrix(inverseMatrix(A));



        double eps = 1e-6;

        System.out.println("Метод половинного ділення: " + bisectionMethod(1, 2, eps));
        System.out.println("Метод хорд: " + chordMethod(1, 2, eps));
        System.out.println("Метод Ньютона: " + newtonMethod(1.5, eps));
        System.out.println("Метод січних: " + secantMethod(1, 2, eps));
        System.out.println("Ітераційний метод: " + iterationMethod(1.5, eps));



        System.out.println("");
        cramerMethod();
        gaussMethod();
        gaussJordanMethod();
        simpleIterationMethod();
        gaussSeidelMethod();



        double a = 0.0;
        double b = Math.PI;
        int n = 10;

        System.out.println("Інтегрування функції f(x) = x^2 * sin(x) на [" + a + ", " + b + "]\n");

        System.out.printf("Метод прямокутників: %.6f\n", rectangleMethod(a, b, n));
        System.out.printf("Метод трапецій: %.6f\n", trapezoidalMethod(a, b, n));
        System.out.printf("Метод Симпсона: %.6f\n", simpsonMethod(a, b, n));
        System.out.printf("Формула Чебишева: %.6f\n", chebyshevMethod(a, b, n));
        System.out.printf("Формула Гауса (2-точкова): %.6f\n", gaussMethod(a, b));



        double x0 = 0.0;
        double y0 = 1.0;
        double h = 0.1;
        int steps = 10;
        System.out.println("");

        System.out.println("Прямий метод Ейлера:");
        printSolution(x0, h, euler(x0, y0, h, steps));

        System.out.println("Зворотний метод Ейлера:");
        printSolution(x0, h, backwardEuler(x0, y0, h, steps));

        System.out.println("Модифікований метод Ейлера:");
        printSolution(x0, h, modifiedEuler(x0, y0, h, steps));

        System.out.println("Виправлений метод Ейлера:");
        printSolution(x0, h, improvedEuler(x0, y0, h, steps));

        System.out.println("Метод Рунге-Кутта:");
        printSolution(x0, h, rungeKutta(x0, y0, h, steps));
    }



    private static double[][] deepCopy(double[][] matrix) {
        double[][] copy = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++)
            System.arraycopy(matrix[i], 0, copy[i], 0, matrix[0].length);
        return copy;
    }



    // а) Складання двох матриць
    public static double[][] addMatrices(double[][] A, double[][] B) {
        int rows = A.length, cols = A[0].length;
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                result[i][j] = A[i][j] + B[i][j];
        return result;
    }

    // а) Віднімання двох матриць
    public static double[][] subtractMatrices(double[][] A, double[][] B) {
        int rows = A.length, cols = A[0].length;
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                result[i][j] = A[i][j] - B[i][j];
        return result;
    }

    // б) Множення матриці на коефіцієнт
    public static double[][] multiplyMatrixByScalar(double[][] A, double scalar) {
        int rows = A.length, cols = A[0].length;
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                result[i][j] = A[i][j] * scalar;
        return result;
    }

    // в) Множення двох матриць
    public static double[][] multiplyMatrices(double[][] A, double[][] B) {
        int rows = A.length;
        int cols = B[0].length;
        int common = B.length;
        double[][] result = new double[rows][cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                for (int k = 0; k < common; k++)
                    result[i][j] += A[i][k] * B[k][j];
        return result;
    }

    // г) Транспонування матриці
    public static double[][] transposeMatrix(double[][] A) {
        int rows = A.length, cols = A[0].length;
        double[][] result = new double[cols][rows];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                result[j][i] = A[i][j];
        return result;
    }

    // д) Знаходження оберненої матриці (для 2x2)
    public static double[][] inverseMatrix(double[][] A) {
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if (det == 0) throw new IllegalArgumentException("Матриця вироджена, не має оберненої.");
        double[][] result = new double[2][2];
        result[0][0] = A[1][1] / det;
        result[0][1] = -A[0][1] / det;
        result[1][0] = -A[1][0] / det;
        result[1][1] = A[0][0] / det;
        return result;
    }

    // Метод для виведення матриці
    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double val : row)
                System.out.printf("%8.2f", val);
            System.out.println();
        }
        System.out.println();
    }



    // Метод половинного ділення
    public static double bisectionMethod(double a, double b, double eps) {
        double c = a;
        while ((b - a) >= eps) {
            c = (a + b) / 2;
            double fc = Math.pow(c, 3) - c - 2;
            double fa = Math.pow(a, 3) - a - 2;
            if (fc == 0.0) break;
            else if (fc * fa < 0) b = c;
            else a = c;
        }
        return c;
    }

    // Метод хорд
    public static double chordMethod(double a, double b, double eps) {
        double x = a;
        while (Math.abs(Math.pow(x, 3) - x - 2) > eps) {
            double fa = Math.pow(a, 3) - a - 2;
            double fb = Math.pow(b, 3) - b - 2;
            x = a - fa * (b - a) / (fb - fa);
            double fx = Math.pow(x, 3) - x - 2;
            if (fx * fa < 0) b = x;
            else a = x;
        }
        return x;
    }

    // Метод Ньютона (дотичних)
    public static double newtonMethod(double x0, double eps) {
        double x1;
        do {
            double fx = Math.pow(x0, 3) - x0 - 2;
            double dfx = 3 * Math.pow(x0, 2) - 1;
            x1 = x0 - fx / dfx;
            if (Math.abs(x1 - x0) < eps) break;
            x0 = x1;
        } while (true);
        return x1;
    }

    // Метод січних
    public static double secantMethod(double x0, double x1, double eps) {
        double x2;
        do {
            double fx0 = Math.pow(x0, 3) - x0 - 2;
            double fx1 = Math.pow(x1, 3) - x1 - 2;
            x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
            if (Math.abs(Math.pow(x2, 3) - x2 - 2) < eps) break;
            x0 = x1;
            x1 = x2;
        } while (true);
        return x2;
    }

    // Ітераційний метод (фіксована ітерація, x = φ(x))
    public static double iterationMethod(double x0, double eps) {
        double x1;
        do {
            // φ(x) = (x + 2)^(1/3) — отримано з рівняння x = (x + 2)^(1/3)
            x1 = Math.cbrt(x0 + 2);
            if (Math.abs(x1 - x0) < eps) break;
            x0 = x1;
        } while (true);
        return x1;
    }



    // Метод Крамера
    public static void cramerMethod() {
        System.out.println("Метод Крамера:");
        double[][] A = {
                {2, -1, 0, 3},
                {1, 0, 2, -1},
                {3, 2, -4, 1},
                {1, -1, 1, -1}
        };
        double[] b = {5, 3, -2, 0};

        double detA = determinant(A);

        if (Math.abs(detA) < 1e-10) {
            System.out.println("Система не має єдиного розв’язку (детермінант = 0)");
            return;
        }

        for (int i = 0; i < 4; i++) {
            double[][] Ai = replaceColumn(A, b, i);
            double detAi = determinant(Ai);
            System.out.printf("x%d = %.4f\n", i + 1, detAi / detA);
        }
    }

    // Метод Гауса
    public static void gaussMethod() {
        System.out.println("Метод Гауса:");
        double[][] A = {
                {2, -1, 0, 3, 5},
                {1, 0, 2, -1, 3},
                {3, 2, -4, 1, -2},
                {1, -1, 1, -1, 0}
        };

        int n = 4;
        for (int i = 0; i < n; i++) {
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }
            double[] temp = A[i];
            A[i] = A[maxRow];
            A[maxRow] = temp;

            for (int k = i + 1; k < n; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j <= n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = A[i][n] / A[i][i];
            for (int k = i - 1; k >= 0; k--) {
                A[k][n] -= A[k][i] * x[i];
            }
        }

        for (int i = 0; i < n; i++) {
            System.out.printf("x%d = %.4f\n", i + 1, x[i]);
        }
    }

    // Метод Гауса-Жордана
    public static void gaussJordanMethod() {
        System.out.println("Метод Гауса-Жордана:");
        double[][] A = {
                {2, -1, 0, 3, 5},
                {1, 0, 2, -1, 3},
                {3, 2, -4, 1, -2},
                {1, -1, 1, -1, 0}
        };

        int n = 4;
        for (int i = 0; i < n; i++) {
            double diag = A[i][i];
            for (int j = 0; j <= n; j++) {
                A[i][j] /= diag;
            }

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = A[k][i];
                    for (int j = 0; j <= n; j++) {
                        A[k][j] -= factor * A[i][j];
                    }
                }
            }
        }

        for (int i = 0; i < n; i++) {
            System.out.printf("x%d = %.4f\n", i + 1, A[i][n]);
        }
    }

    // Ітераційний метод (простих ітерацій)
    public static void simpleIterationMethod() {
        System.out.println("Ітераційний метод (простих ітерацій):");
        int n = 4;
        double[][] A = {
                {10, -1, 2, 0},
                {-1, 11, -1, 3},
                {2, -1, 10, -1},
                {0, 3, -1, 8}
        };
        double[] b = {6, 25, -11, 15};
        double[] x = new double[n];
        double eps = 1e-6;

        boolean converged;
        do {
            converged = true;
            double[] xNew = new double[n];
            for (int i = 0; i < n; i++) {
                double sum = b[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) sum -= A[i][j] * x[j];
                }
                xNew[i] = sum / A[i][i];
                if (Math.abs(xNew[i] - x[i]) > eps) converged = false;
            }
            x = xNew;
        } while (!converged);

        for (int i = 0; i < n; i++) {
            System.out.printf("x%d = %.4f\n", i + 1, x[i]);
        }
    }

    // Ітераційний метод Гауса-Зейделя
    public static void gaussSeidelMethod() {
        System.out.println("Ітераційний метод Гауса-Зейделя:");
        int n = 4;
        double[][] A = {
                {10, -1, 2, 0},
                {-1, 11, -1, 3},
                {2, -1, 10, -1},
                {0, 3, -1, 8}
        };
        double[] b = {6, 25, -11, 15};
        double[] x = new double[n];
        double eps = 1e-6;

        boolean converged;
        do {
            converged = true;
            for (int i = 0; i < n; i++) {
                double old = x[i];
                double sum = b[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) sum -= A[i][j] * x[j];
                }
                x[i] = sum / A[i][i];
                if (Math.abs(x[i] - old) > eps) converged = false;
            }
        } while (!converged);

        for (int i = 0; i < n; i++) {
            System.out.printf("x%d = %.4f\n", i + 1, x[i]);
        }
    }

    private static double determinant(double[][] matrix) {
        int n = matrix.length;

        if (n == 1) {
            return matrix[0][0];
        }

        if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }

        double det = 0.0;

        for (int col = 0; col < n; col++) {
            double[][] subMatrix = new double[n - 1][n - 1];

            for (int i = 1; i < n; i++) {
                int subCol = 0;
                for (int j = 0; j < n; j++) {
                    if (j == col) continue;
                    subMatrix[i - 1][subCol] = matrix[i][j];
                    subCol++;
                }
            }

            det += Math.pow(-1, col) * matrix[0][col] * determinant(subMatrix);
        }

        return det;
    }

    private static double[][] replaceColumn(double[][] A, double[] b, int col) {
        double[][] res = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                res[i][j] = (j == col) ? b[i] : A[i][j];
            }
        }
        return res;
    }



    // Функція для інтегрування: f(x) = x^2 * sin(x)
    private static double f(double x) {
        return x * x * Math.sin(x);
    }

    // Метод прямокутників (лівих)
    // ∫ f(x) dx ≈ h * Σ f(x_i)
    private static double rectangleMethod(double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += f(a + i * h); // лівий прямокутник
        }
        return h * sum;
    }

    // Метод трапецій
    // ∫ f(x) dx ≈ h/2 * [f(x0) + 2f(x1) + ... + 2f(x_{n-1}) + f(xn)]
    private static double trapezoidalMethod(double a, double b, int n) {
        double h = (b - a) / n;
        double sum = (f(a) + f(b)) / 2.0;
        for (int i = 1; i < n; i++) {
            sum += f(a + i * h);
        }
        return h * sum;
    }

    // Метод Симпсона (параболічне наближення)
    // ∫ f(x) dx ≈ h/3 * [f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + f(xn)]
    private static double simpsonMethod(double a, double b, int n) {
        if (n % 2 != 0) n++; // n має бути парним
        double h = (b - a) / n;
        double sum = f(a) + f(b);
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            sum += f(x) * (i % 2 == 0 ? 2 : 4);
        }
        return h / 3.0 * sum;
    }

    // Формула Чебишева (середнє зважене значення на Chebyshev nodes)
    private static double chebyshevMethod(double a, double b, int n) {
        double sum = 0.0;
        for (int k = 1; k <= n; k++) {
            double x = 0.5 * ((b - a) * Math.cos((2.0 * k - 1) / (2.0 * n) * Math.PI) + (b + a));
            sum += f(x);
        }
        return (b - a) * Math.PI / (2.0 * n) * sum;
    }

    // Формула Гауса (2-точкова формула Лежандра)
    // ∫ f(x) dx ≈ (b-a)/2 * [w1*f(x1) + w2*f(x2)]
    private static double gaussMethod(double a, double b) {
        // Вузли та ваги для 2-точкової формули
        double[] x = {-1.0 / Math.sqrt(3), 1.0 / Math.sqrt(3)};
        double[] w = {1.0, 1.0};

        double sum = 0.0;
        for (int i = 0; i < 2; i++) {
            double xi = ((b - a) / 2.0) * x[i] + (a + b) / 2.0;
            sum += w[i] * f(xi);
        }
        return (b - a) / 2.0 * sum;
    }



    // Функція f(x, y) = x + y
    private static double f(double x, double y) {
        return x + y;
    }

    // Прямий метод Ейлера
    public static double[] euler(double x0, double y0, double h, int steps) {
        double[] y = new double[steps + 1];
        y[0] = y0;
        for (int i = 0; i < steps; i++) {
            y[i + 1] = y[i] + h * f(x0 + i * h, y[i]);
        }
        return y;
    }

    // Зворотний метод Ейлера (імпліцитний, чисельне наближення)
    public static double[] backwardEuler(double x0, double y0, double h, int steps) {
        double[] y = new double[steps + 1];
        y[0] = y0;
        for (int i = 0; i < steps; i++) {
            double xNext = x0 + (i + 1) * h;
            double yPrev = y[i];
            // апроксимація методом простого підбору
            double yNext = yPrev;
            for (int j = 0; j < 10; j++) {
                yNext = yPrev + h * f(xNext, yNext);
            }
            y[i + 1] = yNext;
        }
        return y;
    }

    // Модифікований метод Ейлера (середнє значення)
    public static double[] modifiedEuler(double x0, double y0, double h, int steps) {
        double[] y = new double[steps + 1];
        y[0] = y0;
        for (int i = 0; i < steps; i++) {
            double xi = x0 + i * h;
            double k1 = f(xi, y[i]);
            double k2 = f(xi + h, y[i] + h * k1);
            y[i + 1] = y[i] + h * 0.5 * (k1 + k2);
        }
        return y;
    }

    // Виправлений метод Ейлера (ітеративне уточнення)
    public static double[] improvedEuler(double x0, double y0, double h, int steps) {
        double[] y = new double[steps + 1];
        y[0] = y0;
        for (int i = 0; i < steps; i++) {
            double xi = x0 + i * h;
            double predictor = y[i] + h * f(xi, y[i]);
            double corrector = y[i] + h * 0.5 * (f(xi, y[i]) + f(xi + h, predictor));
            y[i + 1] = corrector;
        }
        return y;
    }

    // Метод Рунге-Кутта 4-го порядку
    public static double[] rungeKutta(double x0, double y0, double h, int steps) {
        double[] y = new double[steps + 1];
        y[0] = y0;
        for (int i = 0; i < steps; i++) {
            double xi = x0 + i * h;
            double k1 = h * f(xi, y[i]);
            double k2 = h * f(xi + h / 2, y[i] + k1 / 2);
            double k3 = h * f(xi + h / 2, y[i] + k2 / 2);
            double k4 = h * f(xi + h, y[i] + k3);
            y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        }
        return y;
    }

    // Виведення результатів
    public static void printSolution(double x0, double h, double[] y) {
        for (int i = 0; i < y.length; i++) {
            double x = x0 + i * h;
            System.out.printf("x = %.2f, y = %.6f\n", x, y[i]);
        }
        System.out.println();
    }

}
