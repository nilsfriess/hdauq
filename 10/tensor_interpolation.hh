#pragma once

#include "multiindex.hh"

#include <iostream>
#include <memory>
#include <vector>

template <class OneDInterpolation, class Function,
          std::size_t dim = Function::dim>
struct TensorInterpolation {
  using MI = MultiIndex<dim>;

  TensorInterpolation(const MI &level, const Function &f)
      : oneDOperators(dim), f(f) {

    MI nls;

    for (std::size_t i = 0; i < dim; ++i) {
      oneDOperators[i] = std::make_unique<OneDInterpolation>(level[i]);
      nls[i] = oneDOperators[i]->points.size() - 1;
    }

    indices = computeAllMultiIndices<dim>(nls);
  }

  std::pair<double, std::size_t> evaluate(const std::array<double, dim> &y) {
    double res = 0.0;

    for (const auto &index : indices) {
      std::array<double, dim> evalPoint{0.};
      double basisWeight = 1.;
      for (std::size_t i = 0; i < dim; ++i) {
        auto ind = index[i];

        basisWeight *= oneDOperators[i]->basisFunctions[ind]->operator()(y[i]);
        evalPoint[i] = oneDOperators[i]->points[ind];
      }
      res += basisWeight * f(evalPoint);
    }

    return {res, indices.size()};
  }

  std::vector<std::unique_ptr<OneDInterpolation>> oneDOperators;
  Function f;
  std::vector<MI> indices;
};
