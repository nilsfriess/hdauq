#include "fem.hh"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <cmath>
#include <ctime>
#include <numeric>
#include <vector>

double gen_normal(void) {
  static boost::variate_generator<boost::mt19937, boost::normal_distribution<>>
      generator(boost::mt19937(time(0)), boost::normal_distribution<>());

  double r = generator();
  return r;
}

int main(int /*argc*/, char *argv[]) {
  const double pi = boost::math::constants::pi<double>();

  const int N = 1e3; // Number of MC samples

  const int s = 8;                                 // Where to truncate kappa
  const int n = boost::lexical_cast<int>(argv[1]); // Number of FE mesh points
  const int evalPoint = std::floor(0.7 * n);
  const double q = 2;

  const auto a = [&](int j, double x) {
    return 1 / std::pow(j, q) * std::sin(2 * j * pi * x);
  };

  const auto f = [&](double) { return 1; }; // RHS of BVP

  std::vector<double> results;
  auto bvp = EllipticBVP(f, n);

  std::vector<double> ts(N * s);
  for (int i = 0; i < N * s; ++i) {
    ts[i] = gen_normal();
  }

  for (int k = 0; k < N; ++k) {

    const auto kappa = [&](double x) {
      double res = 0;
      for (int j = 1; j <= s; ++j) {
        res += ts[k * s + j - 1] * a(j, x);
      }

      return std::exp(res);
    };

    bvp.assembleStiffnessMatrix(kappa);
    results.push_back(bvp.solve()[evalPoint]);
  }

  // double sum = std::accumulate(results.begin(), results.end(), 0.0);
  // double mean = sum / results.size();

  for (auto &r : results)
    std::cout << r << "\n";
}
