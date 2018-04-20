#include <iostream>
#include <cmath>
#include "../../euler/src/euler/euler.h"

int main() {
  for (auto i = 0u; i < 3; ++i) {
    /*
     * 7.1
     */
    const auto dt = 1E-3 * std::pow(10, i);
    const auto uo = 1;
    const auto simple = [](const double& t, const double& u) -> double {
      (void) u;
      return -std::sin(t);
    };
    const auto T = 2;

    std::cout << std::endl;
    std::cout << "USING STEP SIZE: " << dt << std::endl;
    std::cout << "==================" << std::endl << std::endl;

    auto soln = Euler::genImplicitEulerSolution<double>(simple, dt, uo);

    std::cout << "7.1:" << std::endl;
    std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
    std::cout << "Exact(2): " << std::cos(T) << std::endl;
    std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;
    std::cout << std::endl;

    /*
     * 7.2
     */
    const auto modified = [](const double& t, const double& u) -> double {
      return -10 * (u - std::cos(t)) - std::sin(t);
    };

    soln = Euler::genImplicitEulerSolution<double>(modified, dt, uo);

    std::cout << "7.2:" << std::endl;
    std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
    std::cout << "Exact(2): " << std::cos(T) << std::endl;
    std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;
    std::cout << std::endl;

    /*
     * 7.3
     */
    const auto largeLambda = [](const double& t, const double& u) -> double {
      return -2100 * (u - std::cos(t)) - std::sin(t);
    };

    soln = Euler::genImplicitEulerSolution<double>(largeLambda, dt, uo);

    std::cout << "7.3:" << std::endl;
    std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
    std::cout << "Exact(2): " << std::cos(T) << std::endl;
    std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;
  }

  return EXIT_SUCCESS;
}
