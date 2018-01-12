#ifndef ABSOLUTE_ERROR_H
#define ABSOLUTE_ERROR_H

#include <cmath>

template <typename T>
/**
 * @param approximate The value that is approximated
 * @param exact The exact value to be compared against
 * @return The relative error
 */
T absoluteError(const T approximate, const T exact) {
  // Absolute error is just the absolute value of the difference
  return std::abs(approximate - exact);
}

#endif
