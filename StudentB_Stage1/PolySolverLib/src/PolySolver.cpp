#include "PolySolver.h"
#include <cmath>
#include <limits>

namespace poly {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

bool validPoly(const double* coeffs, std::size_t n) {
    if (!coeffs || n == 0) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(coeffs[i])) return false;
    }
    return true;
}
} // namespace

double evaluate(const double* coeffs, std::size_t n, double x) {
    if (!validPoly(coeffs, n) || !std::isfinite(x)) return kNaN;
    double acc = coeffs[n - 1];
    for (std::size_t i = n - 1; i > 0; --i) {
        acc = acc * x + coeffs[i - 1];
    }
    return acc;
}

double derivative(const double* coeffs, std::size_t n, double x) {
    if (!validPoly(coeffs, n) || !std::isfinite(x)) return kNaN;
    if (n <= 1) return 0.0;

    double acc = static_cast<double>(n - 1) * coeffs[n - 1];
    for (std::size_t i = n - 1; i > 1; --i) {
        double dk = static_cast<double>(i - 1) * coeffs[i - 1];
        acc = acc * x + dk;
    }
    return acc;
}

int newtonRoot(const double* coeffs, std::size_t n, double x0,
               double tol, int max_iter, double* out_root) {
    if (out_root) *out_root = kNaN;
    if (!validPoly(coeffs, n) || !out_root || !std::isfinite(x0) ||
        !std::isfinite(tol) || tol <= 0.0 || max_iter <= 0) {
        return -1;
    }

    double x = x0;
    for (int iter = 0; iter < max_iter; ++iter) {
        double fx = evaluate(coeffs, n, x);
        if (!std::isfinite(fx)) return -2;
        if (std::fabs(fx) <= tol) {
            *out_root = x;
            return 0;
        }
        double dfx = derivative(coeffs, n, x);
        if (!std::isfinite(dfx) || std::fabs(dfx) < 1e-14) {
            return -2;
        }
        double x_next = x - fx / dfx;
        if (!std::isfinite(x_next)) return -2;
        if (std::fabs(x_next - x) <= tol) {
            *out_root = x_next;
            return 0;
        }
        x = x_next;
    }
    return -2;
}

int findRoots(const double* coeffs, std::size_t n,
              double x_min, double x_max, double step,
              double tol, int max_iter,
              double* out_roots, std::size_t out_cap, std::size_t* out_count) {
    if (out_count) *out_count = 0;
    if (!validPoly(coeffs, n) || !out_roots || !out_count ||
        !std::isfinite(x_min) || !std::isfinite(x_max) || !std::isfinite(step) ||
        step <= 0.0 || x_max < x_min || tol <= 0.0 || max_iter <= 0) {
        return -1;
    }

    std::size_t found = 0;
    for (double cur = x_min; cur <= x_max; cur += step) {
        double root = 0.0;
        if (newtonRoot(coeffs, n, cur, tol, max_iter, &root) == 0) {
            if (std::fabs(evaluate(coeffs, n, root)) <= 1e-5) {
                bool dup = false;
                for (std::size_t k = 0; k < found; ++k) {
                    if (std::fabs(out_roots[k] - root) < tol * 10.0) {
                        dup = true;
                        break;
                    }
                }
                if (!dup) {
                    if (found >= out_cap) {
                        *out_count = found;
                        return -2;
                    }
                    out_roots[found++] = root;
                }
            }
        }
    }
    *out_count = found;
    return 0;
}

} // namespace poly