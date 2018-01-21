#ifndef LOGISTIC_H
#define LOGISTIC_H

#include <functional>
#include <math.h>

template <typename T>
using endomorphism = std::function<T(T)>;

template <typename T>
/**
 * @param alpha Alpha constant in logistic equation
 * @param beta Beta constant in logistic equation
 * @param Po Initial value of equation
 * @return Logistic equation solution
 */
endomorphism<T> genLogistic(T alpha, T beta, T Po) {
  auto amp = ( alpha / Po - beta);

  return endomorphism<T>(
    [=] (T t) {
      return alpha / (amp * exp(-alpha * t) + beta);
    }
  );
}

#endif
