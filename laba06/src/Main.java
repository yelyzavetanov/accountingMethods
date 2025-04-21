public class Main {

    static final double EPSILON = 1e-6;

    public static void main(String[] args) {
        double x = 2.0;
        double y = 1.0;
        int maxIter = 100;

        for (int i = 0; i < maxIter; i++) {
            double f1 = f1(x, y);
            double f2 = f2(x, y);

            double df1dx = df1dx(x, y);
            double df1dy = df1dy(x, y);
            double df2dx = df2dx(x, y);
            double df2dy = df2dy(x, y);

            double det = df1dx * df2dy - df1dy * df2dx;

            if (Math.abs(det) < 1e-10) {
                System.out.println("Матриця Якобі вироджена!");
                return;
            }

            // Метод Крамера
            double dx = (-f1 * df2dy + f2 * df1dy) / det;
            double dy = (-df1dx * f2 + df2dx * f1) / det;

            x += dx;
            y += dy;

            if (Math.abs(dx) < EPSILON && Math.abs(dy) < EPSILON) {
                break;
            }
        }

        System.out.printf("Розв'язок: x = %.6f, y = %.6f%n", x, y);
        System.out.printf("Перевірка: f1 = %.6f, f2 = %.6f%n", f1(x, y), f2(x, y));
    }

    // 🔹 Задаємо саму систему рівнянь
    static double f1(double x, double y) {
        return x * x + y * y - 4;      // x² + y² - 4 = 0
    }

    static double f2(double x, double y) {
        return x * y - 1;              // x * y - 1 = 0
    }

    // 🔹 Похідні
    static double df1dx(double x, double y) {
        return 2 * x;
    }

    static double df1dy(double x, double y) {
        return 2 * y;
    }

    static double df2dx(double x, double y) {
        return y;
    }

    static double df2dy(double x, double y) {
        return x;
    }
}
