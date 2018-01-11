#include "epsilon/epsilon.h"
#include <limits>
#include <iostream>

int main() {
  std::cout << computeEpsilon<double>() << std::endl;
  std::cout << std::numeric_limits<double>::epsilon() << std::endl;

  return EXIT_SUCCESS;
}
