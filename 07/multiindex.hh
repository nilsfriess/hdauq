#pragma once

#include <array>
#include <vector>

template <int dim> using MultiIndex = std::array<unsigned int, dim>;

template <int dim>
std::vector<MultiIndex<dim>> computeAllMultiIndices(const MultiIndex<dim> &m) {
  static_assert(dim > 1, "Error: Dimension must be greater than 1");

  using MultiIndex = MultiIndex<dim>;

  std::vector<MultiIndex> res;

  /* Create all Multiindices that have zero in all but the first directions;
   * in the first direction, we have all values from 0,...,m[0]. */
  for (std::size_t i = 0; i <= m[0]; ++i) {
    MultiIndex mi{0};
    mi[0] = i;
    res.push_back(mi);
  }

  for (std::size_t currDim = 1; currDim < dim; ++currDim) {
    std::vector<MultiIndex> newMultiIndices;
    for (const auto &mi : res) {
      std::size_t currDimCounter = 1; // Case 0 is handeled separately above
      MultiIndex next_mi = mi;
      while (currDimCounter <= m[currDim]) {
        next_mi[currDim] = currDimCounter;
        newMultiIndices.push_back(next_mi);
        ++currDimCounter;
      }
    }
    res.insert(res.end(), newMultiIndices.begin(), newMultiIndices.end());
  }

  return res;
}

template <int dim>
std::vector<MultiIndex<dim>> computeAllMultiIndices(std::size_t l) {
  static_assert(dim > 1, "Error: Dimension must be greater than 1");
  using MultiIndex = MultiIndex<dim>;

  std::vector<MultiIndex> res;

  /* Handle first case separately. This creates all multiindices such that the
   * first two components sum up to l and all other components are zero.
   */
  for (std::size_t i = 0; i <= l; ++i) {
    MultiIndex m{0};
    m[1] = i;
    m[0] = l - i;
    res.push_back(m);
  }

  for (std::size_t currDim = 1; currDim < dim - 1; ++currDim) {
    std::vector<MultiIndex> newMultiIndices;

    for (std::size_t diff = 1; diff <= l; ++diff) {
      for (const auto &mi : res) {
        if (mi[currDim] >= diff) {
          auto next_mi = mi;
          next_mi[currDim] -= diff;
          next_mi[currDim + 1] += diff;
          newMultiIndices.push_back(next_mi);
        }
      }
    }
    res.insert(res.end(), newMultiIndices.begin(), newMultiIndices.end());
  }

  return res;
}
