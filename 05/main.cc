#include <iostream>

#include <boost/lexical_cast.hpp>

#include "fem.hh"

int main(int /*argc*/, char *argv[]) {
  const int n = boost::lexical_cast<int>(argv[1]);

  constexpr auto f = [](double) { return 1; };
  constexpr auto k = [](double) { return 1; };

  EllipticBVP pde(f, k, n);

  // std::cout << pde.getStiffnessMatrix() << std::endl;
  // std::cout << pde.getLoadVector() << std::endl;

  std::cout << pde.solve() << std::endl;
}
