#ifndef RELATIVE_ERROR_H
#define RELATIVE_ERROR_H

#include "../absolute_error/absolute_error.h"

template <typename T>
/**
 * @param approximate The value that is approximated
 * @param exact The exact value to be compared against
 * @return The relative error
 */
T relativeError(const T approximate, const T exact) {
  return absoluteError<T>(approximate, exact) / exact;
}

#endif
