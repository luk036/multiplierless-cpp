#include <cmath>
#include <complex>
#include <vector>

std::vector<double> spectral_fact(std::vector<double> r) {
    // length of the impulse response sequence
    int nr = r.size();
    int n = (nr + 1) / 2;
    // over-sampling factor
    int mult_factor = 30;  // should have mult_factor*(n) >> n
    int m = mult_factor * n;
    // computation method:
    // H(exp(jTw)) = alpha(w) + j*phi(w)
    // where alpha(w) = 1/2*ln(R(w)) and phi(w) = Hilbert_trans(alpha(w))
    // compute 1/2*ln(R(w))
    std::vector<double> w(m);
    for (int i = 0; i < m; i++) {
        w[i] = 2 * M_PI * i / m;
    }
    std::vector<double> R(m);
    for (int i = 0; i < m; i++) {
        for (int j = -(n - 1); j <= n - 1; j++) {
            R[i] += r[j + n - 1] * std::exp(-1i * w[i] * j);
        }
    }
    std::vector<double> alpha(m);
    for (int i = 0; i < m; i++) {
        alpha[i] = 0.5 * std::log(std::abs(R[i]));
    }
    // find the Hilbert transform
    std::vector<std::complex<double>> alphatmp(m);
    for (int i = 0; i < m; i++) {
        alphatmp[i] = std::complex<double>(alpha[i], 0);
    }
    for (int i = floor(m / 2) + 1; i < m; i++) {
        alphatmp[i] = -alphatmp[i];
    }
    alphatmp[0] = 0;
    alphatmp[floor(m / 2) + 1] = 0;
    std::vector<double> phi(m);
    std::vector<std::complex<double>> ifft_input(m);
    for (int i = 0; i < m; i++) {
        ifft_input[i] = std::complex<double>(0, alphatmp[i].imag());
    }
    std::vector<std::complex<double>> ifft_output(m);
    std::vector<double> index;
    for (int i = 0; i < m; i++) {
        if (i % mult_factor == 0) {
            index.push_back(i);
        }
    }
    std::vector<double> alpha1(index.size());
    std::vector<double> phi1(index.size());
    for (int i = 0; i < index.size(); i++) {
        alpha1[i] = alpha[index[i]];
        phi1[i] = phi[index[i]];
    }
    for (int i = 0; i < n; i++) {
        ifft_input[i] = std::exp(std::complex<double>(alpha1[i], phi1[i]));
    }
    for (int i = n + 1; i < m; i++) {
        ifft_input[i] = std::exp(std::complex<double>(alpha1[m - i], -phi1[m - i]));
    }
    for (int i = 0; i < m; i++) {
        ifft_output[i] = ifft_input[i];
    }
    std::vector<double> h(n);
    for (int i = 0; i < n; i++) {
        h[i] = ifft_output[i].real();
    }
    return h;
}
