---
math: true
permalink: /explicitEulerDeux
title: Explicit Euler Part Deux
layout: page
---

**Routine Name**: Testing Explicit Euler

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This code tests the Explicit Euler method implemented in last assignment against the test cases defined in problem 1.

**Input**:

The constants and base cases for the test cases.

**Output**:

This code prints a comparison of the approximated solution (via Explicit Euler) and the exact solution (computed analytically).

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include <array>
#include "../../euler/src/euler/euler.h"
#include "../../testCases/src/testCases/testCases.h"

int main() {
  // Delta t used in Explicit Euler method
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
    const auto approx = Euler::genExplicitEulerSolution<double>(
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
    const auto approx = Euler::genExplicitEulerSolution<double>(
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
approx(0.4) = 14.9182
exact(0.6) = 18.2212
approx(0.6) = 18.2211
exact(0.8) = 22.2554
approx(0.8) = 22.2553
exact(1) = 27.1828
approx(1) = 27.1824

=============
Solving with lambda = -1
exact(0) = 10
approx(0) = 10
exact(0.2) = 8.18731
approx(0.2) = 8.1873
exact(0.4) = 6.7032
approx(0.4) = 6.70319
exact(0.6) = 5.48812
approx(0.6) = 5.4881
exact(0.8) = 4.49329
approx(0.8) = 4.49327
exact(1) = 3.67879
approx(1) = 3.67881

=============
Solving with lambda = 100
exact(0) = 10
approx(0) = 10
exact(0.2) = 4.85165e+09
approx(0.2) = 4.80341e+09
exact(0.4) = 2.35385e+18
approx(0.4) = 2.30727e+18
exact(0.6) = 1.14201e+27
approx(0.6) = 1.10828e+27
exact(0.8) = 5.54062e+35
approx(0.8) = 5.32351e+35
exact(1) = 2.68812e+44
approx(1) = 2.55455e+44

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
approx(0.2) = 22569.9
exact(0.4) = 15815.2
approx(0.4) = 15815
exact(0.6) = 12228
approx(0.6) = 12227.8
exact(0.8) = 10003.8
approx(0.8) = 10003.7
exact(1) = 8490.15
approx(1) = 8490.11
{% endhighlight %}

**Implementation/Code:**

The only new implementation is a comparison function that prints out dialog pertaining to the difference between the exact and found solution.

{% highlight c++ %}
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
{% endhighlight %}
