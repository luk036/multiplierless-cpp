#include <iostream>
#include <cmath>
#include <vector>

std::vector<double> spectral_fact(std::vector<double> r) {
    int n = r.size();
    int mult_factor = 100;
    int m = mult_factor * n;
    std::vector<double> w(m);
    for (int i = 0; i < m; i++) {
        w[i] = 2 * M_PI * i / m;
    }
    std::vector<std::vector<double>> Bn(m, std::vector<double>(n - 1));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n - 1; j++) {
            Bn[i][j] = w[i] * (j + 1);
        }
    }
    std::vector<std::vector<double>> An(m, std::vector<double>(n - 1));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n - 1; j++) {
            An[i][j] = 2 * cos(Bn[i][j]);
        }
    }
    std::vector<std::vector<double>> R(m, std::vector<double>(1));
    for (int i = 0; i < m; i++) {
        R[i][0] = r[0];
        for (int j = 0; j < n - 1; j++) {
            R[i][0] += An[i][j] * r[j + 1];
        }
    }
    std::vector<double> alpha(m);
    for (int i = 0; i < m; i++) {
        alpha[i] = 0.5 * log(abs(R[i][0]));
    }
    std::vector<std::complex<double>> alphatmp(m);
    for (int i = 0; i < m; i++) {
        alphatmp[i] = std::complex<double>(alpha[i], 0);
    }
    int ind = m / 2;
    for (int i = ind; i < m; i++) {
        alphatmp[i] = -alphatmp[i];
    }
    alphatmp[0] = 0;
    alphatmp[ind] = 0;
    std::vector<double> phi(m);
    std::vector<std::complex<double>> tmp(m);
    for (int i = 0; i < m; i++) {
        tmp[i] = 1i * alphatmp[i];
    }
    for (int i = 0; i < m; i++) {
        tmp[i] = exp(tmp[i]);
    }
    for (int i = 0; i < m; i++) {
        tmp[i] *= std::complex<double>(cos(w[i]), sin(w[i]));
    }
    for (int i = 0; i < m; i++) {
        phi[i] = tmp[i].real();
    }
    std::vector<double> alpha1(n);
    std::vector<double> phi1(n);
    for (int i = 0; i < n; i++) {
        alpha1[i] = alpha[i * mult_factor];
        phi1[i] = phi[i * mult_factor];
    }
    std::vector<double> h(n);
    std::vector<std::complex<double>> tmp2(n);
    for (int i = 0; i < n; i++) {
        tmp2[i] = exp(std::complex<double>(alpha1[i], phi1[i]));
    }
    for (int i = 0; i < n; i++) {
        tmp2[i] = std::conj(tmp2[i]);
    }
    for (int i = 0; i < n; i++) {
        tmp2[i] *= h[i];
    }
    for (int i = 0; i < n; i++) {
        h[i] = tmp2[i].real() / n;
    }
    return h;
}

std::vector<double> inverse_spectral_fact(std::vector<double> h) {
    int n = h.size();
    std::vector<double> r(n);
    for (int t = 0; t < n; t++) {
        r[t] = 0;
        for (int i = t; i < n; i++) {
            r[t] += h[i] * h[i - t];
        }
    }
    return r;
}

int main() {
    std::vector<double> r(20);
    for (int i = 0; i < 20; i++) {
        r[i] = (double) rand() / RAND_MAX;
    }
    std::vector<double> h = spectral_fact(r);
    std::vector<double> r2 = inverse_spectral_fact(h);
    return 0;
}
