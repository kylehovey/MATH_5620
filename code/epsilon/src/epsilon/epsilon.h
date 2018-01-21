#ifndef EPSILON_H
#define EPSILON_H

template <typename T>
/**
 * @param eps Initial value for epsilon
 * @return Machine epsilon fir given type
 */
T computeEpsilon(T eps = 1) {
  // While epsilon is still visible to addition
  while (1 + eps != 1) {
    // Reduce exponent by one
    eps /= 2;
  }

  // Eps is no longer visible, show it's value before it was halved
  return eps * 2;
}

#endif
