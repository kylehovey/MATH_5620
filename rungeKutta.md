---
math: true
permalink: /rungeKutta
title: Runge Kutta Methods
layout: page
---

**Routine Name**: Runge Kutta Methods

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This code tests the Runge Kutta Methods of second and fourth order as defined in the book.

**Input**:

The constants and base cases for the test cases.

**Output**:

This code prints a comparison of the approximated solution (via Runge Kutta) and the exact solution (computed analytically).

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include <array>
#include "rungeKutta/rungeKutta.h"
#include "../../testCases/src/testCases/testCases.h"

int main() {
  // Delta t used in Runge Kutta method
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

  std::cout << "|||||||||| Lambda DiffEQ (second order) |||||||||";
  std::cout << std::endl;
  for (const auto lambda : lambdas) {
    const auto approx = RungeKutta::genOrderTwoSolution<double>(
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

  std::cout << "|||||||||| Lambda DiffEQ (fourth order) |||||||||";
  std::cout << std::endl;
  for (const auto lambda : lambdas) {
    const auto approx = RungeKutta::genOrderFourSolution<double>(
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
  std::cout << "|||||||||| Logistic DiffEQ (second order) |||||||||";
  std::cout << std::endl;
  for (const auto Po : Pos) {
    const auto approx = RungeKutta::genOrderTwoSolution<double>(
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

  std::cout << std::endl;
  std::cout << "|||||||||| Logistic DiffEQ (fourth order) |||||||||";
  std::cout << std::endl;
  for (const auto Po : Pos) {
    const auto approx = RungeKutta::genOrderFourSolution<double>(
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
`|||||||||| Lambda DiffEQ (second order) |||||||||

=============
Solving with lambda = 1
exact(0) = 10
approx(0) = 10
exact(0.2) = 12.214
approx(0.2) = 12.214
exact(0.4) = 14.9182
approx(0.4) = 14.9182
exact(0.6) = 18.2212
approx(0.6) = 18.2212
exact(0.8) = 22.2554
approx(0.8) = 22.2554
exact(1) = 27.1828
approx(1) = 27.1825

=============
Solving with lambda = -1
exact(0) = 10
approx(0) = 10
exact(0.2) = 8.18731
approx(0.2) = 8.18731
exact(0.4) = 6.7032
approx(0.4) = 6.7032
exact(0.6) = 5.48812
approx(0.6) = 5.48812
exact(0.8) = 4.49329
approx(0.8) = 4.49329
exact(1) = 3.67879
approx(1) = 3.67883

=============
Solving with lambda = 100
exact(0) = 10
approx(0) = 10
exact(0.2) = 4.85165e+09
approx(0.2) = 4.85164e+09
exact(0.4) = 2.35385e+18
approx(0.4) = 2.35384e+18
exact(0.6) = 1.14201e+27
approx(0.6) = 1.142e+27
exact(0.8) = 5.54062e+35
approx(0.8) = 5.54055e+35
exact(1) = 2.68812e+44
approx(1) = 2.68539e+44
|||||||||| Lambda DiffEQ (fourth order) |||||||||

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
approx(0.2) = 4.9004e+09
exact(0.4) = 2.35385e+18
approx(0.4) = 2.40139e+18
exact(0.6) = 1.14201e+27
approx(0.6) = 1.17677e+27
exact(0.8) = 5.54062e+35
approx(0.8) = 5.76666e+35
exact(1) = 2.68812e+44
approx(1) = 2.82307e+44

|||||||||| Logistic DiffEQ (second order) |||||||||

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
approx(0.2) = 22570.2
exact(0.4) = 15815.2
approx(0.4) = 15815.2
exact(0.6) = 12228
approx(0.6) = 12228
exact(0.8) = 10003.8
approx(0.8) = 10003.8
exact(1) = 8490.15
approx(1) = 8490.22

|||||||||| Logistic DiffEQ (fourth order) |||||||||

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
approx(1) = 8490.32
{% endhighlight %}

**Implementation/Code:**

The Runge Kutta method uses trial steps along the way between each iteration to cancel out lower-order errors. This results in a method that typically converges must faster than other techniques.

{% highlight c++ %}
#include <functional>
#include <cmath>
#include <vector>

namespace RungeKutta {
  template <typename T>
  using endo = std::function<T(const T&)>;

  template <typename T>
  using driver = std::function<T(const T&, const T&)>;

  /**
   * Solve a very basic differential equation using a Runge Kutta technique.
   *  u' = f(t, u)
   * @param {f} RHS of differential equation (utilizes uPrime)
   * @param {dt} Differential of time between samples
   *  (output function will round to these)
   * @param {uInit} Initial value of u(0)
   * @return Function that gives you the output at time t
   */
  template <typename T>
  endo<T> genOrderTwoSolution(
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

          const auto kOne = dt * f(t, lastVal);
          const auto kTwo = dt * f(t + dt / 2, lastVal + kOne / 2);

          cache.push_back(lastVal + kTwo);
        }
      }

      return cache[step];
    };
  }

  /**
   * Solve a very basic differential equation using a Runge Kutta technique.
   *  u' = f(t, u)
   * @param {f} RHS of differential equation (utilizes uPrime)
   * @param {dt} Differential of time between samples
   *  (output function will round to these)
   * @param {uInit} Initial value of u(0)
   * @return Function that gives you the output at time t
   */
  template <typename T>
  endo<T> genOrderFourSolution(
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

          const auto kOne = dt * f(t, lastVal);
          const auto kTwo = dt * f(t + dt / 2, lastVal + kOne / 2);
          const auto kThree = dt * f(t + dt / 2, lastVal + kTwo / 2);
          const auto kFour = dt * f(t + dt, lastVal + kThree);

          cache.push_back(lastVal + kFour);
        }
      }

      return cache[step];
    };
  }
};
{% endhighlight %}
