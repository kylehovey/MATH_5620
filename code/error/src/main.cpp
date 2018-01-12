#include "absolute_error/absolute_error.h"
#include "relative_error/relative_error.h"
#include <limits>
#include <iostream>

int main() {
  const double approx = 3.2;
  const double exact = M_PI;

  std::cout << "PI: " << exact << '\n';
  std::cout << "Approximation: " << approx << '\n';
  std::cout << "- - - - - - - - - -" << '\n';
  std::cout << "Absolute: " << absoluteError<double>(approx, exact) << '\n';
  std::cout << "Relative: " << relativeError<double>(approx, exact) << '\n';

  return EXIT_SUCCESS;
}
