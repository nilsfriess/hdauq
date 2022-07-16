#include <array>
#include <iostream>
#include <random>
#include <vector>

#include "unit_ball_volume.hh"

template <unsigned int dim> int f(const std::array<double, dim> &x) {
  double norm2 = 0.;
  for (unsigned int i = 0; i < dim; ++i)
    norm2 += x[i] * x[i];

  return norm2 <= 1 ? 1 : 0;
}

template <unsigned int dim> double testMC(unsigned int N) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> dist(-1., 1.);
  std::array<double, dim> vec;

  unsigned long sum = 0;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < dim; ++j)
      vec[j] = dist(gen);

    sum += f<dim>(vec);
  }

  double piEstimate = (1 << dim) * double(sum) / N;

  return std::abs(piEstimate - unitBallVolume<dim>());

  // std::cout << "MC Estimate = " << piEstimate << "\n";
  // std::cout << "Exact value = " << unitBallVolume<dim>() << "\n";
  // std::cout << "Error = " << std::abs(piEstimate - unitBallVolume<dim>())
  //           << "\n\n";
}

template <unsigned int dim> void testMCLoop() {
  const unsigned int nMin = 5;
  const unsigned int nMax = 20;

  std::vector<double> results(nMax - nMin);

  unsigned int nRuns = 10;
  for (unsigned int i = 0; i < nRuns; ++i)
    for (unsigned int j = nMin; j < nMax; ++j)
      results[j - nMin] += testMC<dim>(2 << j) / nRuns;

  for (auto res : results)
    std::cout << res << " ";
  std::cout << std::endl;
}

template <unsigned int dim> void printSamples(unsigned int n = 10000) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> dist(-1., 1.);
  std::array<double, dim> vec;

  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < dim; ++j)
      vec[j] = dist(gen);

    std::cout << f<dim>(vec) << " ";
  }

  std::cout << std::endl;
}

template <unsigned int dim> void printQuotientOfVolumes() {
  std::cout << unitBallVolume<dim>() / (1 << dim) << " ";
  printQuotientOfVolumes<dim - 1>();
}

template <> void printQuotientOfVolumes<1>() {}

int main() {
  // (b)
  // testMCLoop<2>();
  // testMCLoop<4>();
  // testMCLoop<6>();
  // testMCLoop<8>();
  // testMCLoop<12>();
  // testMCLoop<14>();
  // testMCLoop<16>();
  // testMCLoop<20>();

  // (c)
  // printSamples<2>();
  // printSamples<3>();
  // printSamples<4>();
  // printSamples<5>();
  // printSamples<6>();
  // printSamples<7>();
  // printSamples<8>();
  // printSamples<9>();
  // printSamples<10>();
  // printSamples<11>();
  // printSamples<12>();

  printQuotientOfVolumes<30>();
}
