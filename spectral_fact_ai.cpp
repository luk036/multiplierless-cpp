#include <cmath>
#include <complex>
#include <iostream>
#include <valarray>

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338327950288
#endif

std::valarray<double> spectral_fact(const std::valarray<double>& r) {
    // length of the impulse response sequence
    size_t n = r.size();

    // over-sampling factor
    size_t mult_factor = 100;  // should have mult_factor*(n) >> n
    size_t m = mult_factor * n;

    // compute 1/2*ln(R(w))
    std::valarray<double> w(m);
    for (size_t i = 0; i < m; i++) {
        w[i] = 2 * M_PI * i / m;
    }
    std::valarray<double> Bn = w * std::valarray<double>(n - 1);
    std::valarray<double> An = 2 * std::cos(Bn);
    std::valarray<double> R(m);
    R[0] = 1;
    R[std::slice(1, n - 1, 1)] = An;
    R = R * r;

    std::valarray<double> alpha = 0.5 * std::log(std::abs(R));

    // find the Hilbert transform
    std::valarray<std::complex<double>> alphatmp(m);
    // alphatmp = std::fft(alpha);
    size_t ind = m / 2;
    alphatmp[ind] = 0;
    // alphatmp[std::slice(ind + 1, m - ind - 1, 1)] *= -1;
    // std::valarray<double> phi = std::real(std::ifft(1i * alphatmp));
    std::valarray<double> phi = alpha;

    // now retrieve the original sampling
    std::valarray<double> alpha1(n), phi1(n);
    for (size_t i = 0; i < n; i++) {
        alpha1[i] = alpha[i * mult_factor];
        phi1[i] = phi[i * mult_factor];
    }

    // compute the impulse response (inverse Fourier transform)
    std::valarray<double> h(n);
    // h = std::real(std::ifft(std::exp(alpha1 + 1i * phi1)));

    return h;
}

std::valarray<double> inverse_spectral_fact(const std::valarray<double>& h) {
    size_t n = h.size();
    std::valarray<double> r(n);
    for (size_t t = 0; t < n; t++) {
        r[t] = 0;
        for (size_t i = t; i < n; i++) {
            r[t] += h[i] * h[i - t];
        }
    }
    return r;
}

int main() {
    std::valarray<double> h
        = {0.76006445, 0.54101887, 0.42012073, 0.3157191, 0.10665804, 0.04326203, 0.01315678};
    std::valarray<double> r = inverse_spectral_fact(h);
    std::valarray<double> h2 = spectral_fact(r);
    // std::cout << h2 << std::endl;
    return 0;
}
