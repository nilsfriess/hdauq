#include <iostream>

#include <array>
#include <numeric>
#include <random>

template <unsigned int n> std::array<double, n> sampleFromX() {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> B(0., 1.);
  std::exponential_distribution T(1.);

  std::array<double, n> res;
  for (unsigned int i = 0; i < n; ++i)
    res[i] = B(gen) / T(gen);

  return res;
}

template <unsigned int n> double computeMCEstimate() {
  auto samples = sampleFromX<n>();
  auto res = std::accumulate(samples.begin(), samples.end(), 0.0);
  return res / n;
}

int main() {
  for (int i = 0; i < 10; ++i) {
    std::cout << computeMCEstimate<2 << 5>() << " ";
    std::cout << computeMCEstimate<2 << 6>() << " ";
    std::cout << computeMCEstimate<2 << 7>() << " ";
    std::cout << computeMCEstimate<2 << 8>() << " ";
    std::cout << computeMCEstimate<2 << 9>() << " ";
    std::cout << computeMCEstimate<2 << 10>() << " ";
    std::cout << computeMCEstimate<2 << 11>() << " ";
    std::cout << computeMCEstimate<2 << 12>() << " ";
    std::cout << computeMCEstimate<2 << 13>() << " ";
    std::cout << computeMCEstimate<2 << 14>() << " ";
    std::cout << computeMCEstimate<2 << 15>() << " ";
    std::cout << std::endl;
  }
}
