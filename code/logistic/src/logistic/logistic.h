#ifndef LOGISTIC_H
#define LOGISTIC_H

#include <functional>
#include <math.h>

template <typename T>
std::function<T(T)> genLogistic(T alpha, T beta, T Po) {
  auto amp = ( alpha / Po - beta);

  return std::function<T(T)>([amp, alpha, beta](T t) {
    return alpha / (amp * exp(-alpha * t) + beta);
  });
}

#endif
