#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

std::vector<double> spectral_fact(std::vector<double> r) {
    size_t n = r.size();
    size_t mult_factor = 100;
    size_t m = mult_factor * n;
    std::vector<double> w(m);
    for (size_t i = 0; i < m; i++) {
        w[i] = 2 * M_PI * double(i) / double(m);
    }
    std::vector<std::vector<double>> Bn(m, std::vector<double>(n - 1));
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n - 1; j++) {
            Bn[i][j] = w[i] * double(j + 1);
        }
    }
    std::vector<std::vector<double>> An(m, std::vector<double>(n - 1));
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n - 1; j++) {
            An[i][j] = 2 * cos(Bn[i][j]);
        }
    }
    std::vector<std::vector<double>> R(m, std::vector<double>(1));
    for (size_t i = 0; i < m; i++) {
        R[i][0] = r[0];
        for (size_t j = 0; j < n - 1; j++) {
            R[i][0] += An[i][j] * r[j + 1];
        }
    }
    std::vector<double> alpha(m);
    for (size_t i = 0; i < m; i++) {
        alpha[i] = 0.5 * log(abs(R[i][0]));
    }
    std::vector<std::complex<double>> alphatmp(m);
    for (size_t i = 0; i < m; i++) {
        alphatmp[i] = std::complex<double>(alpha[i], 0);
    }
    size_t ind = m / 2;
    for (size_t i = ind; i < m; i++) {
        alphatmp[i] = -alphatmp[i];
    }
    alphatmp[0] = 0;
    alphatmp[ind] = 0;
    std::vector<double> phi(m);
    std::vector<std::complex<double>> tmp(m);
    for (size_t i = 0; i < m; i++) {
        tmp[i] = std::complex<double>(0.0, 1.0) * alphatmp[i];
    }
    for (size_t i = 0; i < m; i++) {
        tmp[i] = exp(tmp[i]);
    }
    for (size_t i = 0; i < m; i++) {
        tmp[i] *= std::complex<double>(cos(w[i]), sin(w[i]));
    }
    for (size_t i = 0; i < m; i++) {
        phi[i] = tmp[i].real();
    }
    std::vector<double> alpha1(n);
    std::vector<double> phi1(n);
    for (size_t i = 0; i < n; i++) {
        alpha1[i] = alpha[i * mult_factor];
        phi1[i] = phi[i * mult_factor];
    }
    std::vector<double> h(n);
    std::vector<std::complex<double>> tmp2(n);
    for (size_t i = 0; i < n; i++) {
        tmp2[i] = exp(std::complex<double>(alpha1[i], phi1[i]));
    }
    for (size_t i = 0; i < n; i++) {
        tmp2[i] = std::conj(tmp2[i]);
    }
    for (size_t i = 0; i < n; i++) {
        tmp2[i] *= h[i];
    }
    for (size_t i = 0; i < n; i++) {
        h[i] = tmp2[i].real() / double(n);
    }
    return h;
}

std::vector<double> inverse_spectral_fact(std::vector<double> h) {
    size_t n = h.size();
    std::vector<double> r(n);
    for (size_t t = 0; t < n; t++) {
        r[t] = 0;
        for (size_t i = t; i < n; i++) {
            r[t] += h[i] * h[i - t];
        }
    }
    return r;
}

int main() {
    std::vector<double> r(20);
    for (size_t i = 0; i < 20; i++) {
        r[i] = (double)rand() / RAND_MAX;
    }
    std::vector<double> h = spectral_fact(r);
    std::vector<double> r2 = inverse_spectral_fact(h);
    return 0;
}
