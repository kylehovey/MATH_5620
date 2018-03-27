#include <iostream>
#include "euler/euler.h"

int main() {
  // Domain and range of solution
  using D = double;

  // Driving function and initial condition
  const auto uInit = 1;
  const auto dt = 0.01;
  const Euler::driver<D> f = [](const D& t, const D& ut) -> D {
    return t * ut;
  };

  // Exact solution
  const Euler::endo<D> exact = [](const D& t) -> D {
    return std::exp(t * t / 2);
  };

  // Solve the system
  const auto soln = Euler::genEulerSolution<D>(f, dt, uInit);

  // Test the solution
  for (auto i = 0; i < 10; ++i) {
    std::cout << "u(" << i * dt << ")" << std::endl;
    std::cout << exact(i * dt) << " vs ";
    std::cout << soln(i * dt) << std::endl;
  }

  return EXIT_SUCCESS;
}
