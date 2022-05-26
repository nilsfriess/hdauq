#include <iostream>

#include <boost/lexical_cast.hpp>

#include "fem.hh"

int main(int /*argc*/, char *argv[]) {
  const int n = boost::lexical_cast<int>(argv[1]);

  constexpr auto f = [](double) { return 1; };
  constexpr auto k = [](double x) { return 1 / (2 + x); };

  EllipticBVP pde(f, k, n);

  // std::cout << pde.getStiffnessMatrix() << std::endl;
  // std::cout << pde.getLoadVector() << std::endl;

  auto res = pde.solve();
  for (auto entry : res)
    std::cout << entry << "\n";
}
