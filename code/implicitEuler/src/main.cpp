#include <iostream>
#include <array>
#include "../../euler/src/euler/euler.h"
#include "../../testCases/src/testCases/testCases.h"

/**
 * Compare the output of two endomorphisms over a given
 * domain
 * @param exact First function
 * @param approx Second function
 * @param range The beginning and ending of the range
 * @param steps Steps to iterate over
 */
template <typename T>
void printComparison(
    const Euler::endo<T>& exact,
    const Euler::endo<T>& approx,
    const std::tuple<T, T> range,
    const unsigned int& steps
) {
  const auto [ start, stop ] = range;
  const auto dt = (stop - start) / steps;

  for (auto i = 0u; i <= steps; ++i) {
    const auto t = start + i * dt;
    std::cout << "exact(" << t << ") = " << exact(t) << std::endl;
    std::cout << "approx(" << t << ") = " << approx(t) << std::endl;
  }
}

int main() {
  // Delta t used in Implicit Euler method
  const auto dt = 0.00001;

  // Variables for evaluation
  const std::tuple<double, double> domain = { 0.0, 1.0 };
  const unsigned int steps = 5;

  // Test cases for u' = λu
  const auto alpha = 10.0;
  const std::array<double, 3> lambdas = { 1, -1, 100 };

  // Test cases for Logistic Equation
  const auto gamma = 0.1;
  const auto beta = 0.0001;
  const std::array<double, 2> Pos = { 25, 40000 };

  std::cout << "|||||||||| Lambda DiffEQ |||||||||" << std::endl;
  for (const auto lambda : lambdas) {
    const auto approx = Euler::genImplicitEulerSolution<double>(
        [=](const double& t, const double& u) -> double {
          (void) t;
          return lambda * u;
        },
        dt,
        alpha
    );

    const auto exact = TestCases::genLambdaSolution<double>(lambda, alpha);

    std::cout << std::endl << "=============" << std::endl;
    std::cout << "Solving with lambda = " << lambda << std::endl;

    printComparison<double>(exact, approx, domain, steps);
  }

  std::cout << std::endl;
  std::cout << "|||||||||| Logistic DiffEQ |||||||||" << std::endl;
  for (const auto Po : Pos) {
    const auto approx = Euler::genImplicitEulerSolution<double>(
        [=](const double& t, const double& P) -> double {
          (void) t;

          return gamma * P - beta * P * P;
        },
        dt,
        Po
    );

    const auto exact = TestCases::genLogisticSolution(beta, gamma, Po);

    std::cout << std::endl << "=============" << std::endl;
    std::cout << "Solving with Po = " << Po << std::endl;

    printComparison<double>(exact, approx, domain, steps);
  }

  return EXIT_SUCCESS;
}
