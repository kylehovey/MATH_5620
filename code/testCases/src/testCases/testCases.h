#ifndef TEST_CASES_H
#define TEST_CASES_H

#include <functional>
#include <cmath>

template <typename T>
using endomorphism = std::function<T(const T&)>;

namespace TestCases {
  /**
   * Generates the analytic solution to the differential equation
   *  u' = λu; u(0) := α
   * @param lambda Lambda in the equation above
   * @param alpha Initial condition at u(0)
   * @return Analytic solution
   */
  template <typename T>
  endomorphism<T> genLambdaSolution(const T& lambda, const T& alpha) {
    return [=](const T& t) -> T {
      return alpha * std::exp(lambda * t);
    };
  }

  /**
   * Generates the analytic solution to the logistic equation
   *  P' = γP - βP; P(0) := Po
   * @param gamma Gamma in above equation
   * @param beta Beta in above equation
   * @param Po Initial condition at P(0)
   * @return Analytic solution
   */
  template <typename T>
  endomorphism<T> genLogisticSolution(
      const T& beta,
      const T& gamma,
      const T& Po
  ) {
    const auto amp = ( gamma / Po - beta);

    return [=](const T& t) {
      return gamma / (amp * exp(-gamma * t) + beta);
    };
  }
};

#endif
