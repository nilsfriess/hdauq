#pragma once

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <boost/math/quadrature/trapezoidal.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

using Mat = Eigen::SparseMatrix<double>;
using Vec = Eigen::VectorXd;

using boost::math::quadrature::trapezoidal;

template <typename F, typename K> class EllipticBVP {

public:
  EllipticBVP(F f_, K k_, int n_)
      : f(f_), k(k_), n(n_), x(n + 2), h(1. / (n + 1)) {
    // Assemble Grid
    x[0] = 0;
    for (int i = 1; i < n + 2; ++i)
      x[i] = i * h;
    // std::cout << x[x.size() - 1] << std::endl;
    x[x.size() - 1] = 1; // Ensure last point is 1

    assembleStiffnessMatrix();
    assembleLoadVector();
  }

  Vec getLoadVector() const { return b; }
  Mat getStiffnessMatrix() const { return A; }

  Vec solve() const {
    Eigen::SimplicialLDLT<Mat, Eigen::Lower> solver;
    return solver.compute(A).solve(b);
  }

private:
  /* Assembles the sparse stiffness matrix, but only the
   * lower diagonal and the main diagonal, since it
   * is symmetric.
   */
  void assembleStiffnessMatrix() {

    // Construct lower diagonal
    Mat U{n, n};
    U.reserve(n - 1);

    for (int i = 2; i < n + 1; ++i) {
      auto k_dphi = [&](double z) {
        return k(z) * d_phi(i - 1, z) * d_phi(i, z);
      };

      auto res = trapezoidal(k_dphi, x[i - 1] + 10e-12, x[i] - 10e-12);
      U.insert(i - 1, i - 2) = res;

      // std::cout << "Res = " << res << "Error = " << std::abs(-1 * (n + 1) -
      // res)
      //           << std::endl;

      // U.insert(i - 1, i - 2) = res;
      // if (i < n)
      //   U.insert(i - 1, i - 2) = -1 * (n + 1) + 0.06;
      // else
      // U.insert(i, i - 1) = -1 * (n + 1);
    }

    // std::cout << U;

    // Construct diagonal
    Mat D{n, n};
    D.reserve(n);

    for (int i = 1; i < n + 1; ++i) {
      auto k_dphi = [&](double z) { return k(z) * d_phi(i, z) * d_phi(i, z); };
      D.insert(i - 1, i - 1) = trapezoidal(k_dphi, x[i - 1], x[i + 1]);
    }

    A = (D + U);
  }

  void assembleLoadVector() {
    b = Vec(n);

    for (int i = 1; i < n + 1; ++i) {
      auto f_phi = [&](double z) { return f(z) * phi(i, z); };
      b(i - 1) = trapezoidal(f_phi, x[i - 1], x[i + 1]);
    }
  }

  double phi(int i, double val) {
    if (val < x[i])
      return (val - x[i - 1]) / h;
    else
      return (x[i + 1] - val) / h;
  }

  double d_phi(int i, double val) {
    if (val < x[i])
      return 1. / h;
    else
      return -1. / h;
  }

  F f;
  K k;

  int n;

  std::vector<double> x;
  const double h;

  Mat A; // Stiffness matrix
  Vec b; // Load vector
};
