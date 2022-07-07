#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "hat_function.hh"
#include "one_d_interpolation.hh"
#include "sparse_grid_interpolation.hh"
#include "tensor_interpolation.hh"

template <std::size_t s> struct Function {
  static const std::size_t dim = s;

  double operator()(const std::array<double, dim> &points) {
    double sum = 0.0;

    for (std::size_t i = 0; i < dim; ++i)
      sum += (points[i] - 0.5) / ((i + 1) * (i + 1));

    return 1. / (1. + sum);
  }
};

template <std::size_t dim> void testFunc(std::size_t level) {
  using Func = Function<dim>;
  Func f;

  std::array<double, dim> y;
  for (std::size_t i = 0; i < dim; ++i)
    y[i] = 1 / std::sqrt(2.);

  for (std::size_t l = 2; l <= level; ++l) {
    SparseGridInterpolation<OneDInterpolation<HatFunction>, Func> interpolation(
        l, f);

    auto [res, funcEvals] = interpolation.evaluate(y);

    std::cout << "Level: " << l << "\n";
    std::cout << "Exact value:  " << std::setprecision(12) << f(y) << "\n";
    std::cout << "Interpolated: " << res << "\n";
    std::cout << "Function evaluations: " << funcEvals << "\n\n";
  }
}

int main() {
  constexpr std::size_t dim = 10;

  testFunc<dim>(5);
}
