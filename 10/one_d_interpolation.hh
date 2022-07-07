#pragma once

#include <memory>
#include <vector>

template <class BasisFunction> struct OneDInterpolation {
  explicit OneDInterpolation(std::size_t level)
      : n((1 << level) + 1), points(n), basisFunctions(n) {
    for (std::size_t i = 0; i < n; ++i)
      points[i] = double(i) / (n - 1);

    for (std::size_t i = 0; i < n; ++i) {
      basisFunctions[i] = std::make_unique<BasisFunction>(i, points);
    }
  }

  std::size_t n;
  std::vector<double> points;
  std::vector<std::unique_ptr<BasisFunction>> basisFunctions;
};
