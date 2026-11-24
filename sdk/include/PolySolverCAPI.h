#pragma once
#include <stddef.h>

#if defined(_WIN32) || defined(_WIN64)
#ifdef POLYSOLVER_EXPORTS
#define POLY_API extern "C" __declspec(dllexport)
#else
#define POLY_API extern "C" __declspec(dllimport)
#endif
#else
#define POLY_API extern "C"
#endif

POLY_API double poly_evaluate(const double* coeffs, size_t n, double x);
POLY_API double poly_derivative(const double* coeffs, size_t n, double x);
POLY_API int poly_newton_root(const double* coeffs, size_t n, double x0,
    double tol, int max_iter, double* out_root);
POLY_API int poly_find_roots(const double* coeffs, size_t n,
    double x_min, double x_max, double step,
    double tol, int max_iter,
    double* out_roots, size_t out_cap, size_t* out_count);