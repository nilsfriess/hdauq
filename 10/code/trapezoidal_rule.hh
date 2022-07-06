#pragma once

#include <algorithm>
#include <vector>

struct TrapezoidalRule {
  explicit TrapezoidalRule(std::size_t nl, double a = 0., double b = 1.)
      : nl(nl), points(nl + 1), weights(nl + 1) {
    double h = (b - a) / nl;

    for (std::size_t i = 0; i < nl + 1; ++i)
      points[i] = a + i * h;

    std::fill(weights.begin() + 1, weights.end() - 1, h);
    weights[0] = h / 2;
    weights[nl] = h / 2;
  }

  std::size_t nl;
  std::vector<double> points;
  std::vector<double> weights;
};
