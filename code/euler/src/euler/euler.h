#ifndef EULER_H
#define EULER_H

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
  endo<T> genExplicitEulerSolution(
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

          cache.push_back(lastVal + dt * f(i * dt, lastVal));
        }
      }

      return cache[step];
    };
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
                return (nextU - lastVal) / dt - f(i * dt, nextU);
              },
              lastVal
          ));
        }
      }

      return cache[step];
    };
  }
};

#endif
