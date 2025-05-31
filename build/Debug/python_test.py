import sys
import os

# Upewnij się, że folder z .pyd znajduje się w sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, "build", "Debug"))  # lub gdziekolwiek jest .pyd


import example  # Twój moduł C++ z pybind11

print("=== Generator sygnałów + DFT ===")
print("Wybierz typ sygnału do narysowania:")
print("sin    - Sygnał sinusoidalny")
print("cos    - Sygnał cosinusoidalny")
print("square - Sygnał prostokątny")
print("saw    - Sygnał piłokształtny")

signal_type = input("Wpisz typ sygnału: ").strip().lower()

# Parametry
frequency = float(input("Podaj częstotliwość (Hz): "))
amplitude = float(input("Podaj amplitudę: "))
phase = float(input("Podaj fazę (w radianach): "))
samples = int(input("Podaj liczbę próbek: "))

# Generuj sygnał i weź dane
if signal_type == "sin":
    signal = example.generate_sine(frequency, amplitude, phase, samples)
elif signal_type == "cos":
    signal = example.generate_cosine(frequency, amplitude, phase, samples)
elif signal_type == "square":
    signal = example.generate_square(frequency, amplitude, phase, samples)
elif signal_type == "saw":
    signal = example.generate_sawtooth(frequency, amplitude, phase, samples)
else:
    print("Nieznany typ sygnału.")
    exit()

# Rysowanie oryginalnego sygnału
example.plot_signal(signal, "Oryginalny sygnał")

# DFT
dft_result = example.compute_dft(signal)
example.plot_magnitude(dft_result, "Widmo DFT (moduł)")

# IDFT
reconstructed = example.compute_idft(dft_result)
example.plot_signal(reconstructed, "Sygnał po IDFT")



# Wykryj krawędzie z progiem czułości
edges = example.detect_edges(signal, threshold=0.1)

# Rysuj z krawędziami na czerwono
example.plot_signal_with_edges(signal, edges, "Sygnał z wykrytymi krawędziami")



