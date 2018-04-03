---
math: true
permalink: /implicitEuler
title: Implicit Euler
layout: page
---

**Routine Name**: Testing Implicit Euler

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This code tests the Implicit Euler method.

**Input**:

The constants and base cases for the test cases.

**Output**:

This code prints a comparison of the approximated solution (via Implicit Euler) and the exact solution (computed analytically).

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include <array>
#include "../../euler/src/euler/euler.h"
#include "../../testCases/src/testCases/testCases.h"

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
    const auto approx = Euler::genImplicit<double>(
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
    const auto approx = Euler::genImplicit<double>(
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
{% endhighlight %}

Output:

{% highlight C++ %}
|||||||||| Lambda DiffEQ |||||||||

=============
Solving with lambda = 1
exact(0) = 10
approx(0) = 10
exact(0.2) = 12.214
approx(0.2) = 12.214
exact(0.4) = 14.9182
approx(0.4) = 14.9183
exact(0.6) = 18.2212
approx(0.6) = 18.2212
exact(0.8) = 22.2554
approx(0.8) = 22.2555
exact(1) = 27.1828
approx(1) = 27.1827

=============
Solving with lambda = -1
exact(0) = 10
approx(0) = 10
exact(0.2) = 8.18731
approx(0.2) = 8.18732
exact(0.4) = 6.7032
approx(0.4) = 6.70321
exact(0.6) = 5.48812
approx(0.6) = 5.48813
exact(0.8) = 4.49329
approx(0.8) = 4.49331
exact(1) = 3.67879
approx(1) = 3.67885

=============
Solving with lambda = 100
exact(0) = 10
approx(0) = 10
exact(0.2) = 4.85165e+09
approx(0.2) = 4.90044e+09
exact(0.4) = 2.35385e+18
approx(0.4) = inf
exact(0.6) = 1.14201e+27
approx(0.6) = inf
exact(0.8) = 5.54062e+35
approx(0.8) = inf
exact(1) = 2.68812e+44
approx(1) = inf

|||||||||| Logistic DiffEQ |||||||||

=============
Solving with Po = 25
exact(0) = 25
approx(0) = 25
exact(0.2) = 25.4922
approx(0.2) = 25.4922
exact(0.4) = 25.9937
approx(0.4) = 25.9937
exact(0.6) = 26.5049
approx(0.6) = 26.5049
exact(0.8) = 27.0259
approx(0.8) = 27.0259
exact(1) = 27.5568
approx(1) = 27.5568

=============
Solving with Po = 40000
exact(0) = 40000
approx(0) = 40000
exact(0.2) = 22570.2
approx(0.2) = 22570.4
exact(0.4) = 15815.2
approx(0.4) = 15815.4
exact(0.6) = 12228
approx(0.6) = 12228.2
exact(0.8) = 10003.8
approx(0.8) = 10004
exact(1) = 8490.15
approx(1) = 8490.32 highlight C++ %}
{% endhighlight %}

**Implementation/Code:**

A Newton method approach is used to solve implicitly for the next value in the iteration as a function of the previous value. Just as in the Explicit method, this method returns a memoized version of the function that closes over a call cache so that the solution is dynamically created.

{% highlight c++ %}
#include <functional>
#include <cmath>
#include <vector>

namespace Euler {
  template <typename T>
  using endo = std::function<T(const T&)>;

  template <typename T>
  using driver = std::function<T(const T&, const T&)>;

  /**
   * Use newton's method to find a local zero of a function given a guess
   * @param f Function to find the zero of
   * @param guess Initial guess for the location of the zero
   * @param tolerance Absolute error away from zero tolerated
   * @param dt Time step to use for secant approximation of derivative
   * @param maxIterations Maximum iterations before exiting
   * @return An input x such that |f(x)| < tolerance
   */
  template <typename T>
  T newton(
      const endo<T>& f,
      const T& guess,
      const T& tolerance = 0.001,
      const T& dt = 0.001,
      const unsigned int& maxIter = 100
  ) {
    auto out = guess;

    for (auto i = 0u; std::abs(f(out)) > tolerance && i < maxIter; ++i) {
      const auto df = (f(out + dt) - f(out)) / dt;
      out = out - f(out) / df;
    }

    return out;
  }

  /**
   * Solve a very basic differential equation using an Explicit Euler technique.
   *  u' = f(t, u)
   * @param {f} RHS of differential equation (utilizes uPrime)
   * @param {dt} Differential of time between samples
   *  (output function will round to these)
   * @param {uInit} Initial value of u(0)
   * @return Function that gives you the output at time t
   */
  template <typename T>
  endo<T> genImplicitEulerSolution(
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

          //cache.push_back(lastVal + dt * f(t, lastVal));
          cache.push_back(newton<double>(
              [&](const double& nextU) -> double {
                return (nextU - lastVal) / dt - f(t, nextU);
              },
              lastVal
          ));
        }
      }

      return cache[step];
    };
  }
};
{% endhighlight %}
