#ifndef LOGISTIC_H
#define LOGISTIC_H

#include <functional>
#include <complex>
#include <iostream>

namespace SpringMass {
  // We will be using this type a lot, make an alias for it
  template <typename T>
  using endomorphism = std::function<T(T)>;

  // Make an alias for the zero function
  template <typename T>
  const endomorphism<T> zeroFunction = endomorphism<T>(
    /**
     * @param t Temporal input
     * @return Zero cast into the template type
     */
    [](T t) {
      (void) t;
      return (T) 0;
    }
  );

  template <typename T>
  /**
   * Generate an analytic solution to a spring mass system that gives
   * the position of the mass at any time t.
   * @param yo Initial position (t = 0)
   * @param dyo Initial derivative (t = 0)
   * @param m The mass
   * @param gamma Damping constant (air resistance, etc...)
   * @param k Spring constant
   * @param endomorphism f Driving function (non-homogenous part)
   * @return Solution to spring-mass equation
   */
  endomorphism<T> genSpringMass(
      T yo = 0,
      T dyo = 0,
      T m = 1,
      T gamma = 1,
      T k = 1,
      endomorphism<T> f = zeroFunction<T>
  ) {
    // Compute the eigenvalues and allow for complex solutions
    auto desc = std::sqrt(
        (std::complex<T>) (std::pow(gamma, 2) - 4 * m * k)
    );

    std::complex<T> eigens[2] = {
      (-gamma + desc) / std::complex<double>(2 * m),
      (-gamma - desc) / std::complex<double>(2 * m)
    };

    // Determine whether or not we have real or repeated roots
    auto distinct = eigens[0] != eigens[1];

    // Determine the leading coefficients of the solution
    auto denom = eigens[1] - eigens[0];

    std::complex<T> coeffs[2] = {
      distinct ? (eigens[1] * yo - dyo) / denom : yo,
      distinct ? (-eigens[0] * yo + dyo) / denom : dyo - eigens[0] * yo
    };

    // Determine additional u(t) required for linear independance
    const endomorphism<T> u(
      distinct ?
        (endomorphism<T>) [](T t) { (void) t; return 1; } :
        (endomorphism<T>) [](T t) { return t; }
    );

    // Stub for determining particular solution
    // TODO: Find particular solution
    (void) f;
    endomorphism<T> particular = zeroFunction<T>;

    return endomorphism<T>(
      /**
       * @param t Temporal input
       * @return Position of mass in system at t
       */
      [=] (T t) {
        /*
         * Since the characteristic equation has real coefficients, both
         * eigenvalues will be complex conjugates of each other. Since
         * conjugation is preserved under exponentiation, we are adding
         * two conjugates and thus cancelling the imaginary part. This
         * means that we only have to take the real part of the solution
         * (assuming that the particular solution has real range).
         */
        return (
          coeffs[0] * std::exp(eigens[0] * t) +
          coeffs[1] * std::exp(eigens[1] * t) * u(t) +
          particular(t)
        ).real();
      }
    );
  }
}

#endif
