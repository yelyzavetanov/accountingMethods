public class Main {

    // Функція правої частини рівняння y'(x) = x + x^2 + 0.1 * x^3
    public static double f(double x, double y) {
        return x + Math.pow(x, 2) + 0.1 * Math.pow(x, 3);
    }

    // Метод Ейлера
    public static void eulerMethod(double x0, double y0, double h, double xn) {
        double x = x0, y = y0;
        System.out.println("Метод Ейлера:");
        while (x <= xn) {
            System.out.printf("x = %.2f, y = %.5f\n", x, y);
            y = y + h * f(x, y);
            x += h;
        }
    }

    // Метод Рунге-Кутта 4-го порядку
    public static void rungeKuttaMethod(double x0, double y0, double h, double xn) {
        double x = x0, y = y0;
        System.out.println("\nМетод Рунге-Кутта:");
        while (x <= xn) {
            System.out.printf("x = %.2f, y = %.5f\n", x, y);

            double k1 = h * f(x, y);
            double k2 = h * f(x + h / 2, y + k1 / 2);
            double k3 = h * f(x + h / 2, y + k2 / 2);
            double k4 = h * f(x + h, y + k3);

            y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            x += h;
        }
    }

    public static void main(String[] args) {
        double x0 = 0, y0 = 0, h = 0.1, xn = 1;
        eulerMethod(x0, y0, h, xn);
        rungeKuttaMethod(x0, y0, h, xn);
    }
}
