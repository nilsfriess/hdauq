#pragma once

#include <cstddef>
#include <vector>

#include <iostream>

struct HatFunction {
  HatFunction(std::size_t j, const std::vector<double> &points)
      : j(j), points(points) {}

  double operator()(double y) {
    auto nl = points.size() - 1;

    if ((j > 0) and (y >= points[j - 1]) and (y <= points[j]))
      return (y - points[j - 1]) / (points[j] - points[j - 1]);
    else if ((j < nl) and (y >= points[j]) and (y <= points[j + 1]))
      return (points[j + 1] - y) / (points[j + 1] - points[j]);
    else
      return 0.;
  }

  std::size_t j;
  std::vector<double> points;
};
