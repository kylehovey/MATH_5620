#ifndef PREDICTOR_CORRECTOR_H
#define PREDICTOR_CORRECTOR_H

#include <vector>
#include <array>
#include <functional>
#include <cmath>

namespace PredCorr {
  template <typename T>
  using endo = std::function<T(const T&)>;

  template <typename T>
  using driver = std::function<T(const T&, const T&)>;
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
          const auto tNow = i * dt;

          const auto kOne = lastVal + dt * f(tNow, lastVal);
          const auto kTwo = lastVal + 0.5 * dt *
            (f(tNow, lastVal) + f(tNow, kOne));

          cache.push_back(kTwo);
        }
      }

      return cache[step];
    };
  }
};

#endif
