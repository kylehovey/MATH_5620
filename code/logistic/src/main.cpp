#include "logistic/logistic.h"
#include <limits>
#include <iostream>

int main() {
  // Parameters for logistic equation
  const double
    alpha = 2,
    beta = 1,
    Po = 1;

  // Generate the required function
  auto logistic = Logistic::genLogistic<double>(alpha, beta, Po);

  // Call it for some basic values
  for (int i = -10; i <= 10; i++) {
    std::cout << i << " -> " << logistic(i) << std::endl;
  }

  return EXIT_SUCCESS;
}
