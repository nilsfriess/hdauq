#include <cmath>
#include <iostream>

#include "gauss_legendre.hh"
#include "sparse_grid_quadrature.hh"

template <typename T, std::size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &arr) {
  for (const auto &elem : arr)
    out << elem << " ";
  return out;
}

template <int s> struct TestFunc {
  static const int dim = s;

  double operator()(std::array<double, dim> y) const {
    double res = 1;
    double factor = (2 + 1 / (2. * dim));
    for (int i = 0; i < dim; ++i)
      res *= factor;

    for (int i = 0; i < dim; ++i)
      res *= std::pow(y[i], 1 + 1 / (2. * dim));
    return res;
  }
};

template <std::size_t dim> void testProblem(std::size_t l) {
  using Function = TestFunc<dim>;
  SparseGridQuadratureRule<GaussLegendre, Function> quad(l, Function());

  auto [res, funcEvals] = quad.integrate();

  std::cout << dim << " " << l << " " << std::setprecision(12) << res << " "
            << funcEvals << std::endl;
}

int main() {
  testProblem<5>(1);
  testProblem<5>(2);
  testProblem<5>(3);
  testProblem<5>(4);
  testProblem<5>(5);
  testProblem<5>(6);

  testProblem<6>(1);
  testProblem<6>(2);
  testProblem<6>(3);
  testProblem<6>(4);
  testProblem<6>(5);
  testProblem<6>(6);

  testProblem<7>(1);
  testProblem<7>(2);
  testProblem<7>(3);
  testProblem<7>(4);
  testProblem<7>(5);
  testProblem<7>(6);

  testProblem<8>(1);
  testProblem<8>(2);
  testProblem<8>(3);
  testProblem<8>(4);
  testProblem<8>(5);
  testProblem<8>(6);

  // auto indices = computeAllMultiIndices<6>(4);

  // for (const auto &ind : indices) {
  //   std::cout << "( " << ind << ")\n";
  // }
  // std::cout << "Total: " << indices.size() << "\n";
}
