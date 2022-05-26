#pragma once

#include <Eigen/Sparse>

#include <boost/math/quadrature/trapezoidal.hpp>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

#include "thomas.hh"

using boost::math::quadrature::trapezoidal;

template <typename F> class EllipticBVP {

public:
  EllipticBVP(F f_, int n_) : f(f_), n(n_), x(n + 2), h(1. / (n + 1)) {
    // Assemble Grid
    x[0] = 0;
    for (int i = 1; i < n + 2; ++i)
      x[i] = i * h;
    // std::cout << x[x.size() - 1] << std::endl;
    x[x.size() - 1] = 1; // Ensure last point is 1

    assembleLoadVector();
  }

  std::vector<double> getLoadVector() const { return b; }
  std::array<std::vector<double>, 2> getStiffnessMatrix() const { return A; }

  std::vector<double> solve() const {
    return thomas_solve(A[0], A[1], A[0], b);
  }

  /* Assembles the sparse stiffness matrix, but only the
   * lower diagonal and the main diagonal, since it
   * is symmetric.
   */
  template <typename K> void assembleStiffnessMatrix(K k) {

    std::vector<double> lowerDiag(n - 1);
    for (int i = 2; i < n + 1; ++i) {
      auto k_dphi = [&](double z) {
        return k(z) * d_phi(i - 1, z) * d_phi(i, z);
      };

      auto res = trapezoidal(k_dphi, x[i - 1] + 10e-12, x[i] - 10e-12);
      lowerDiag[i - 2] = res;
    }

    std::vector<double> mainDiag(n);

    for (int i = 1; i < n + 1; ++i) {
      auto k_dphi = [&](double z) { return k(z) * d_phi(i, z) * d_phi(i, z); };
      mainDiag[i - 1] = trapezoidal(k_dphi, x[i - 1], x[i + 1]);
    }

    A[0] = lowerDiag;
    A[1] = mainDiag;
  }

  void assembleLoadVector() {
    b = std::vector<double>(n);

    for (int i = 1; i < n + 1; ++i) {
      auto f_phi = [&](double z) { return f(z) * phi(i, z); };
      b[i - 1] = trapezoidal(f_phi, x[i - 1], x[i + 1]);
    }
  }

private:
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

  int n;

  std::vector<double> x;
  const double h;

  std::array<std::vector<double>, 2> A; // Stiffness matrix
  std::vector<double> b;                // Load vector
};
