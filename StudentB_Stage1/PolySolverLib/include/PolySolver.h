#pragma once
#include <cstddef>

namespace poly {

// Значення полінома P(x) = coeffs[0] + coeffs[1]*x + ... + coeffs[n-1]*x^(n-1)
double evaluate(const double* coeffs, std::size_t n, double x);

// Значення похідної P'(x) у точці x
double derivative(const double* coeffs, std::size_t n, double x);

// Пошук одного кореня методом Ньютона зі стартового наближення x0
// Повертає: 0 — корінь знайдено, -1 — некоректні параметри, -2 — не збіглося / похідна ~ 0
int newtonRoot(const double* coeffs, std::size_t n, double x0,
               double tol, int max_iter, double* out_root);

// Пошук дійсних коренів на відрізку [x_min, x_max] із кроком step
// Повертає: 0 — успіх, -1 — некоректні параметри, -2 — замалий буфер out_roots
int findRoots(const double* coeffs, std::size_t n,
              double x_min, double x_max, double step,
              double tol, int max_iter,
              double* out_roots, std::size_t out_cap, std::size_t* out_count);

} // namespace poly