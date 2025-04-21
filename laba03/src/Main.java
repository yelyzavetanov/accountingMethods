public class Main {
    static final double EPSILON = 1e-6; // Точність

    // Функція f(x) = x + sin(x)
    static double f(double x) {
        return x + Math.sin(x);
    }

    // Метод хорд
    static double chordMethod(double a, double b) {
        double x;
        do {
            x = a - f(a) * (b - a) / (f(b) - f(a));
            if (f(a) * f(x) < 0) {
                b = x;
            } else {
                a = x;
            }
        } while (Math.abs(f(x)) > EPSILON);
        return x;
    }

    // Метод січних
    static double secantMethod(double x0, double x1) {
        double x2;
        do {
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
            x0 = x1;
            x1 = x2;
        } while (Math.abs(f(x2)) > EPSILON);
        return x2;
    }

    public static void main(String[] args) {
        double a = -2, b = 0; // Початковий інтервал
        double x0 = -2, x1 = -1; // Початкові наближення для січних

        System.out.println("Корінь методом хорд: " + chordMethod(a, b));
        System.out.println("Корінь методом січних: " + secantMethod(x0, x1));
    }
}
