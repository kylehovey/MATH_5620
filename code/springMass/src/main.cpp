#include "springMass/springMass.h"
#include <limits>
#include <iostream>

int main() {
  // Constants and initial values
  const double
      m = 5,
      gamma = 1,
      k = 10,
      yo = 0,
      dyo = 5;

  // Generate the analytic solution
  auto y = genSpringMass<double>(yo, dyo, m, gamma, k);

  // Output position at time t for 3 seconds (damped oscillation)
  for (double t = 0; t < 3; t += 0.1) {
    std::cout << y(t) << std::endl;
  }

  return EXIT_SUCCESS;
}
