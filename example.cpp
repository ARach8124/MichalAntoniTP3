#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <matplot/matplot.h>
#include <cmath>
#include <complex>
#include <vector>

namespace py = pybind11;
namespace plt = matplot;

// Stała pi
const double M_PI = 3.141592653589793;

// Funkcja generująca wektor x dla sygnałów
std::vector<double> generate_x(double frequency, int samples = 1000) {
    double period = 1.0 / frequency;
    return plt::linspace(0, 2 * period, samples);
}

// Funkcja do rysowania sygnału
void plot_signal(const std::vector<double>& signal, const std::string& title) {
    std::vector<double> x = generate_x(1.0); // Zakładając częstotliwość 1Hz na wykresie
    plt::plot(x, signal);
    plt::title(title);
    plt::xlabel("Czas [s]");
    plt::ylabel("Amplituda");
    plt::show();
}

// Funkcja do rysowania widma DFT
void plot_magnitude(const std::vector<std::complex<double>>& dft_result, const std::string& title) {
    std::vector<double> magnitudes;
    for (const auto& value : dft_result) {
        magnitudes.push_back(std::abs(value)); // Obliczamy moduł DFT
    }

    std::vector<double> x(magnitudes.size());
    for (size_t i = 0; i < magnitudes.size(); ++i) {
        x[i] = i;
    }

    plt::plot(x, magnitudes);
    plt::title(title);
    plt::xlabel("Częstotliwość");
    plt::ylabel("Moduł");
    plt::show();
}

// Funkcja obliczająca DFT
std::vector<std::complex<double>> compute_dft(const std::vector<double>& signal) {
    int N = signal.size();
    std::vector<std::complex<double>> dft(N);

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0, 0);
        for (int n = 0; n < N; ++n) {
            double angle = 2 * M_PI * k * n / N;
            sum += signal[n] * std::exp(std::complex<double>(0, -angle));
        }
        dft[k] = sum;
    }
    return dft;
}

// Funkcja obliczająca odwrotną DFT
std::vector<double> compute_idft(const std::vector<std::complex<double>>& dft_result) {
    int N = dft_result.size();
    std::vector<double> idft(N);

    for (int n = 0; n < N; ++n) {
        std::complex<double> sum(0, 0);
        for (int k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            sum += dft_result[k] * std::exp(std::complex<double>(0, angle));
        }
        idft[n] = sum.real() / N; // Dzielenie przez N, aby uzyskać normalizację
    }

    return idft;
}

// Funkcja generująca sygnał sinusoidalny
std::vector<double> generate_sine(double frequency, double amplitude, double phase, int samples) {
    auto x = generate_x(frequency, samples);
    std::vector<double> y;
    for (double t : x) {
        y.push_back(amplitude * sin(2 * M_PI * frequency * t + phase));
    }
    return y;
}

// Funkcja generująca sygnał cosinusoidalny
std::vector<double> generate_cosine(double frequency, double amplitude, double phase, int samples) {
    auto x = generate_x(frequency, samples);
    std::vector<double> y;
    for (double t : x) {
        y.push_back(amplitude * cos(2 * M_PI * frequency * t + phase));
    }
    return y;
}

// Funkcja generująca sygnał prostokątny
std::vector<double> generate_square(double frequency, double amplitude, double phase, int samples) {
    auto x = generate_x(frequency, samples);
    std::vector<double> y;
    for (double t : x) {
        double value = sin(2 * M_PI * frequency * t + phase);
        y.push_back(value >= 0 ? amplitude : -amplitude);
    }
    return y;
}

// Funkcja generująca sygnał piłokształtny
std::vector<double> generate_sawtooth(double frequency, double amplitude, double phase, int samples) {
    auto x = generate_x(frequency, samples);
    std::vector<double> y;
    for (double t : x) {
        double time_shifted = t + phase / (2 * M_PI * frequency);
        double fractional = time_shifted * frequency - floor(time_shifted * frequency);
        y.push_back(amplitude * (2 * fractional - 1)); // od -A do A
    }
    return y;
}

std::vector<size_t> detect_edges(const std::vector<double>& signal, double threshold) {
    std::vector<size_t> edge_indices;
    for (size_t i = 1; i < signal.size() - 1; ++i) {
        double diff = std::abs(signal[i + 1] - signal[i - 1]);
        if (diff > threshold) {
            edge_indices.push_back(i);
        }
    }
    return edge_indices;
}
void plot_signal_with_edges(const std::vector<double>& signal,
                            const std::vector<size_t>& edge_indices,
                            const std::string& title) {
    std::vector<double> x = generate_x(1.0, signal.size());
    plt::plot(x, signal);  // sygnał

    std::vector<double> edge_x, edge_y;
    for (size_t i : edge_indices) {
        edge_x.push_back(x[i]);
        edge_y.push_back(signal[i]);
    }

    plt::hold(true);
    plt::plot(edge_x, edge_y, "r*");  // czerwone gwiazdki
    plt::title(title);
    plt::xlabel("Czas [s]");
    plt::ylabel("Amplituda");
    plt::show();
}




PYBIND11_MODULE(example, m) {
    m.def("generate_sine", &generate_sine, "Generuje sygnał sinusoidalny", py::arg("frequency"), py::arg("amplitude"), py::arg("phase"), py::arg("samples"));
    m.def("generate_cosine", &generate_cosine, "Generuje sygnał cosinusoidalny", py::arg("frequency"), py::arg("amplitude"), py::arg("phase"), py::arg("samples"));
    m.def("generate_square", &generate_square, "Generuje sygnał prostokątny", py::arg("frequency"), py::arg("amplitude"), py::arg("phase"), py::arg("samples"));
    m.def("generate_sawtooth", &generate_sawtooth, "Generuje sygnał piłokształtny", py::arg("frequency"), py::arg("amplitude"), py::arg("phase"), py::arg("samples"));
    m.def("plot_signal", &plot_signal, "Rysuje sygnał", py::arg("signal"), py::arg("title"));
    m.def("plot_magnitude", &plot_magnitude, "Rysuje widmo DFT", py::arg("dft_result"), py::arg("title"));
    m.def("compute_dft", &compute_dft, "Oblicza DFT", py::arg("signal"));
    m.def("compute_idft", &compute_idft, "Oblicza IDFT", py::arg("dft_result"));
    m.def("detect_edges", &detect_edges, "Wykrywa indeksy krawędzi", py::arg("signal"), py::arg("threshold"));
    m.def("plot_signal_with_edges", &plot_signal_with_edges, "Rysuje sygnał z zaznaczonymi krawędziami", py::arg("signal"), py::arg("edge_indices"), py::arg("title"));


}