#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using ComplexNum = std::complex<double>;

// Function to calculate the Euclidean norm (magnitude) of a complex number
double magnitude(const ComplexNum& c) {
    return std::sqrt(c.real() * c.real() + c.imag() * c.imag());
}

std::vector<ComplexNum> lejaOrder(std::vector<ComplexNum> points) {
    // Check if input is empty
    if (points.empty()) {
        throw std::runtime_error("Input must be a non-empty vector of points.");
    }

    // Sort points by magnitude initially
    std::sort(points.begin(), points.end(),
              [](const ComplexNum& a, const ComplexNum& b) { return magnitude(a) < magnitude(b); });

    std::vector<ComplexNum> lejaOrderedPoints;
    lejaOrderedPoints.push_back(points.front());  // Start with the smallest magnitude point
    points.erase(points.begin());                 // Remove this point from further consideration

    while (!points.empty()) {
        // Find the point with the maximum minimum distance to the current Leja sequence
        size_t nextIdx = 0;
        double maxMinDist = 0.0;
        for (size_t i = 0; i < points.size(); ++i) {
            double minDist = std::numeric_limits<double>::max();
            for (const auto& p : lejaOrderedPoints) {
                minDist = std::min(minDist, magnitude(points[i] - p));
            }
            if (minDist > maxMinDist) {
                maxMinDist = minDist;
                nextIdx = i;
            }
        }

        // Append this point to the Leja ordered sequence
        lejaOrderedPoints.push_back(points[nextIdx]);
        // Remove this point from further consideration
        points.erase(points.begin() + nextIdx);
    }

    return lejaOrderedPoints;
}

int main() {
    std::vector<ComplexNum> points = {{1.0, 1.0}, {2.0, -2.0}, {0.5, 0.5}, {-1.0, 0.0}, {3.0, 3.0}};
    std::vector<ComplexNum> lejaOrdered = lejaOrder(points);

    for (const auto& point : lejaOrdered) {
        std::cout << point << std::endl;
    }

    return 0;
}
