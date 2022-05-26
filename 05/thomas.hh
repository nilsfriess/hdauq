#pragma once

#include <vector>

using Vec = std::vector<double>;

Vec thomas_solve(Vec lower, Vec main, Vec upper, Vec rhs) {
  int n = main.size() - 1;

  lower.insert(lower.begin(), 0);

  upper[0] /= main[0];
  rhs[0] /= main[0];

  for (int i = 1; i < n; i++) {
    upper[i] /= main[i] - lower[i] * upper[i - 1];
    rhs[i] =
        (rhs[i] - lower[i] * rhs[i - 1]) / (main[i] - lower[i] * upper[i - 1]);
  }

  rhs[n] =
      (rhs[n] - lower[n] * rhs[n - 1]) / (main[n] - lower[n] * upper[n - 1]);

  for (int i = n; i-- > 0;)
    rhs[i] -= upper[i] * rhs[i + 1];

  return rhs;
}
