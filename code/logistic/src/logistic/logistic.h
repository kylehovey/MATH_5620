#ifndef LOGISTIC_H
#define LOGISTIC_H

#include <functional>
#include <math.h>

template <typename T>
using homomorphism = std::function<T(T)>;

template <typename T>
homomorphism<T> genLogistic(T alpha, T beta, T Po) {
  auto amp = ( alpha / Po - beta);

  return homomorphism<T>(
    [=] (T t) {
      return alpha / (amp * exp(-alpha * t) + beta);
    }
  );
}

#endif
