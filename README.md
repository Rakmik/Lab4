# Контракт бібліотеки PolySolverDLL (Варіант 6: Розв'язання поліноміальних рівнянь)

- **Версія API:** 1.0
- **Угода зв'язування:** `extern "C"` (C-ABI).
- **Формат полінома:** `coeffs[0] + coeffs[1]*x + ... + coeffs[n-1]*x^(n-1)` довжини `n`.
- **Керування пам'яттю:** буфер `out_roots` виділяє та звільняє викликаюча сторона[cite: 24].

## Експортовані функції

1. `double poly_evaluate(const double* coeffs, size_t n, double x)`
   - **Призначення:** обчислює значення полінома в точці `x` за схемою Горнера.
   - **Повертає:** значення `P(x)` або `NaN` при некоректних аргументах.

2. `double poly_derivative(const double* coeffs, size_t n, double x)`
   - **Призначення:** обчислює першу похідну `P'(x)` в точці `x`.
   - **Повертає:** `0.0` для константи (`n == 1`), `P'(x)` для вищих степенів або `NaN` при помилці.

3. `int poly_newton_root(const double* coeffs, size_t n, double x0, double tol, int max_iter, double* out_root)`
   - **Призначення:** шукає корінь методом Ньютона зі стартової точки `x0` із точністю `tol` та лімітом `max_iter`.
   - **Повертає:** `0` — корінь записано в `*out_root`; `-1` — некоректні параметри (`coeffs == nullptr`, `n == 0`, `out_root == nullptr`, `tol <= 0`, `max_iter <= 0`, `NaN`/`Inf`); `-2` — не збіглося або `|P'(x)| < 1e-14`.

4. `int poly_find_roots(const double* coeffs, size_t n, double x_min, double x_max, double step, double tol, int max_iter, double* out_roots, size_t out_cap, size_t* out_count)`
   - **Призначення:** сканує `[x_min, x_max]` з кроком `step` та уточнює корені методом Ньютона.
   - **Повертає:** `0` — успіх; `-1` — некоректні аргументи (`step <= 0`, `x_max < x_min`, `nullptr` тощо); `-2` — недостатня місткість `out_cap`. Записує кількість знайдених коренів у `*out_count`.
