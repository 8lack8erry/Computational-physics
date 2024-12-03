from matplotlib import pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq  # Import fft and fftfreq from scipy.fft

d = np.loadtxt("osAsf.txt")
t_po = d[:,0]
x_po = d[:,1]

# t = t_po
# x = x_po

t = []
x = []
min = 0
max = 50

for i in range(len(t_po)):
    if t_po[i] > min and t_po[i] < max:
        t.append(t_po[i])
        x.append(x_po[i])

fx = fft(x)
sum_fx = sum(np.abs(fx))

sampling_frequency = 1 / (t[1] - t[0])
freq = fftfreq(len(fx), 1 / sampling_frequency)

normalized_amplitude = np.abs(fx) / sum_fx

plt.stem(freq, normalized_amplitude, 'b', markerfmt=" ", basefmt="-b", label=r"$fft, t\in(0, 50)$")
plt.xlabel(r"$\omega$")
plt.ylabel("fft")
plt.xlim(0, 2)
plt.axvline(x=1 / (3 * np.pi), color='r', linestyle='--', label=r"$\omega_{Forzante}$")
plt.legend()
plt.savefig("FFT_ArmonicOscillator.png")