#include <iostream>
#include <vector>
#include "PolySolver.h"

int main() {
    // Поліном P(x) = (x - 1)(x - 2)(x + 3) = -6 + 5x + 2x^2 - x^3
    // Дійсні корені: -3, 1, 2
    std::vector<double> coeffs = { -6.0, 5.0, 2.0, -1.0 };

    std::cout << "P(1)  = " << poly::evaluate(coeffs.data(), coeffs.size(), 1.0) << "\n";
    std::cout << "P'(1) = " << poly::derivative(coeffs.data(), coeffs.size(), 1.0) << "\n";

    double root = 0.0;
    if (poly::newtonRoot(coeffs.data(), coeffs.size(), 1.2, 1e-9, 100, &root) == 0) {
        std::cout << "Newton root near 1.2: " << root << "\n";
    }

    std::vector<double> roots(10, 0.0);
    std::size_t count = 0;
    int rc = poly::findRoots(coeffs.data(), coeffs.size(), -5.0, 5.0, 0.5, 1e-9, 100,
        roots.data(), roots.size(), &count);

    std::cout << "Found roots (rc = " << rc << ", count = " << count << "): ";
    for (std::size_t i = 0; i < count; ++i) std::cout << roots[i] << " ";
    std::cout << "\n";

    return 0;
}