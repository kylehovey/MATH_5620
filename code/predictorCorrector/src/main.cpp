#include <iostream>
#include <vector>
#include <array>
#include "../../testCases/src/testCases/testCases.h"

template <typename T>
using endo = std::function<T(const T&)>;

template <typename T>
using driver = std::function<T(const T&, const T&)>;

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
    const endo<T>& exact,
    const endo<T>& approx,
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

/**
 * Solve a basic differential equation using a Predictor Corrector technique.
 *  u' = f(t, u)
 * @param {f} RHS of differential equation (utilizes uPrime)
 * @param {dt} Differential of time between samples
 *  (output function will round to these)
 * @param {uInit} Initial value of u(0)
 * @return Function that gives you the output at time t
 */
template <typename T>
endo<T> predictorCorrector(
    const driver<T>& f,
    const T& dt,
    const T& uInit
) {
  // Memoization cache
  std::vector<T> cache = { };
  cache.push_back(uInit);

  return [=](const T& t) mutable -> T {
    const auto step = std::floor(t / dt);

    if (step >= cache.size()) {
      const auto size = cache.size();

      for (auto i = size; i <= step; ++i) {
        const auto lastVal = cache[i - 1];

        const auto kOne = lastVal + dt * f(t, lastVal);
        const auto kTwo = lastVal + 0.5 * dt * (f(t, lastVal) + f(t, kOne));

        cache.push_back(kTwo);
      }
    }

    return cache[step];
  };
}

int main() {
  // Delta t used in Predictor Corrector method
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
    const auto approx = predictorCorrector<double>(
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
    const auto approx = predictorCorrector<double>(
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
