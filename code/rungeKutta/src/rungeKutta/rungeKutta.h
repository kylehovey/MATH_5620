#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

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

#endif
