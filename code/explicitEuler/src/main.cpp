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
    std::cout << "u(" << i * dt << ")" << std::endl;
    std::cout << exact(i * dt) << " vs ";
    std::cout << soln(i * dt) << std::endl;
  }

  return EXIT_SUCCESS;
}
