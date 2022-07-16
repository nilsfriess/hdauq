#pragma once

#include <cstdint>

constexpr double pi = 3.14159265358979323846;

template <unsigned int s>
constexpr double unitBallVolume() {
  return 2*pi/s*unitBallVolume<s-2>();
}


template <>
constexpr double unitBallVolume<0>() {
  return 1;
}

template <>
constexpr double unitBallVolume<1>() {
  return 2;
}
