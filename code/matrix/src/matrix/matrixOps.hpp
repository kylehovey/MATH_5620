#ifndef MATRIX_OPS_HPP
#define MATRIX_OPS_HPP

#include "matrix.h"

namespace Matrix {
  template <typename U>
  Matrix<U> operator+(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.add(rhs);
  }

  template <typename U>
  Matrix<U> operator-(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.subtract(rhs);
  }

  template <typename U, typename V>
  Matrix<V> operator*(const U& lhs, const Matrix<V>& rhs) {
    return rhs.scalarMult(lhs);
  }

  template <typename U, typename V>
  Matrix<V> operator*(const Matrix<V>& lhs, const U& rhs) {
    return lhs.scalarMult(rhs);
  }

  template <typename U>
  bool operator==(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.isEqualTo(rhs);
  }

  template <typename U>
  bool operator!=(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return !(lhs == rhs);
  }
}

#endif
