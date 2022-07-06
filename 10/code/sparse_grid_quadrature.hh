#pragma once

#include <boost/math/special_functions/binomial.hpp>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "multiindex.hh"
#include "tensor_quad.hh"

template <class OneDRule, class Function, std::size_t dim = Function::dim>
struct SparseGridQuadratureRule {
  using TensorRule = TensorQuadratureRule<OneDRule, Function>;

  SparseGridQuadratureRule(std::size_t level, const Function &f)
      : level(level), f(f) {}

  std::pair<double, std::size_t> integrate() const {
    double res = 0.0;
    std::size_t funcEvals = 0;

    for (std::size_t r = 0; r <= std::min(level, dim - 1); ++r) {
      double fac = boost::math::binomial_coefficient<double>(dim - 1, r);
      if (r % 2 != 0)
        fac *= -1;

      double inner_res = 0.0;
      auto indices = computeAllMultiIndices<dim>(level - r);
      for (const auto &ind : indices) {
        auto tensorRule = TensorQuadratureRule<OneDRule, Function>(ind, f);
        auto [next_res, next_evals] = tensorRule.integrate();
        inner_res += next_res;
        funcEvals += next_evals;
      }
      res += fac * inner_res;
    }
    return {res, funcEvals};
  }

  std::size_t level;
  Function f;
};
