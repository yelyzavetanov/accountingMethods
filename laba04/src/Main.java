import java.util.function.Function;

public class Main {
    // Функція для інтегрування
    static double function(double x) {
        return (x * x) / Math.log(x);
    }

    // Метод прямокутників (Ньютона-Котеса)
    static double rectangleMethod(Function<Double, Double> f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            double x = a + i * h + h / 2; // Середина підінтервалу
            sum += f.apply(x);
        }
        return sum * h;
    }

    // Метод Гауса для 2 точок
    static double gaussMethod(Function<Double, Double> f, double a, double b) {
        double[] nodes = {-0.5773502692, 0.5773502692}; // Корені полінома Лежандра для n=2
        double[] weights = {1.0, 1.0}; // Вага для n=2

        double sum = 0;
        for (int i = 0; i < nodes.length; i++) {
            double x = 0.5 * (b - a) * nodes[i] + 0.5 * (a + b);
            sum += weights[i] * f.apply(x);
        }
        return 0.5 * (b - a) * sum;
    }

    public static void main(String[] args) {
        double a = 1, b = 10;
        int n = 50;

        // Метод прямокутників
        double resultRectangle = rectangleMethod(Main::function, a, b, n);
        System.out.println("Інтеграл методом прямокутників: " + resultRectangle);

        // Метод Гауса (інтегруємо на кожному з 5 підінтервалів)
        double resultGauss = 0;
        double step = (b - a) / n;
        for (int i = 0; i < n; i++) {
            resultGauss += gaussMethod(Main::function, a + i * step, a + (i + 1) * step);
        }
        System.out.println("Інтеграл методом Гауса: " + resultGauss);
    }
}
