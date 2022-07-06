#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>

#include "multiindex.hh"

template <class OneDRule, class Function, int dim = Function::dim>
struct TensorQuadratureRule {
  using MI = MultiIndex<dim>;

  TensorQuadratureRule(const MI &level, const Function &f)
      : oneDRules(dim), f(f) {
    // Map level vector to vector of number of points
    std::transform(level.cbegin(), level.cend(), nls.begin(),
                   [](std::size_t l) { return 2 << l; });

    for (int i = 0; i < dim; ++i)
      oneDRules[i] = std::make_unique<OneDRule>(nls[i]);

    indices = computeAllMultiIndices<dim>(nls);
  }

  std::pair<double, std::size_t> integrate() {
    double res = 0;
    for (const auto &index : indices) {
      double weight = 1;
      std::array<double, dim> evalPoint{0.};
      for (int i = 0; i < dim; ++i) {
        auto ind = index[i];
        weight *= oneDRules[i]->weights[ind];
        evalPoint[i] = oneDRules[i]->points[ind];
      }

      res += weight * f(evalPoint);
    }
    return {res, indices.size()};
  }

  MI nls;
  std::vector<std::unique_ptr<OneDRule>> oneDRules;
  Function f;
  std::vector<MI> indices;
};
