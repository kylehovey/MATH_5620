---
math: true
permalink: /explicitEuler
title: Explicit Euler
layout: page
---

**Routine Name**: Explicit Euler Method

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This computational method solves an initial value problem of the form

# \\[y'(t) = f(t, y(t)) \wedge y(t_0) \equiv y_0\\]

**Input**:

The driving function, initial condition, and time step.

**Output**:

A function that will evaluate to an approximated solution. The returned function also memoizes entries and uses dynamic programming to never compute the same value twice.

**Usage/Example**:

{% highlight C++ %}
int main() {
  // Domain and range of solution
  using D = double;

  // Driving function and initial condition
  const auto uInit = 1;
  const auto dt = 0.01;
  const driver<D> f = [](const D& t, const D& ut) -> D {
    return t * ut;
  };

  // Exact solution
  const endo<D> exact = [](const D& t) -> D {
    return std::exp(t * t / 2);
  };

  // Solve the system
  const auto soln = genEulerSolution<D>(f, dt, uInit);

  // Test the solution
  for (auto i = 0; i < 10; ++i) {
    std::cout << "u(" << i * dt << ")" << " vs ";
    std::cout << exact(i * dt) << std::endl;
    std::cout << soln(i * dt) << std::endl;
  }

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
u(0)
1 vs 1
u(0.01)
1.00005 vs 1.0001
u(0.02)
1.0002 vs 1.0003
u(0.03)
1.00045 vs 1.0006
u(0.04)
1.0008 vs 1.001
u(0.05)
1.00125 vs 1.0015
u(0.06)
1.0018 vs 1.0021
u(0.07)
1.00245 vs 1.0028
u(0.08)
1.00321 vs 1.00361
u(0.09)
1.00406 vs 1.00451
{% endhighlight %}

**Implementation/Code:**

{% highlight c++ %}
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

template <typename T>
using endo = std::function<T(const T&)>;

template <typename T>
using driver = std::function<T(const T&, const T&)>;

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
endo<T> genEulerSolution(
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

        cache.push_back(lastVal + dt * f(t, lastVal));
      }
    }

    return cache[step];
  };
}
{% endhighlight %}
