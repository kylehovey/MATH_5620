#ifndef MATRIX_OPS_HPP
#define MATRIX_OPS_HPP

#include "matrix.h"
#include <iostream>

namespace Matrix {
  template <typename U>
  Matrix<U> operator+(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.add(rhs);
  }

  template <typename U>
  Matrix<U> operator-(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.subtract(rhs);
  }

  template <typename U>
  Matrix<U> operator*(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.multiply(rhs);
  }

  template <typename U, typename V>
  Matrix<U> operator*(const V& lhs, const Matrix<U>& rhs) {
    return rhs.scalarMult((U) lhs);
  }

  template <typename U, typename V>
  Matrix<U> operator*(const Matrix<U>& lhs, const V& rhs) {
    return lhs.scalarMult((U) rhs);
  }

  template <typename U>
  bool operator==(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return lhs.isEqualTo(rhs);
  }

  template <typename U>
  bool operator!=(const Matrix<U>& lhs, const Matrix<U>& rhs) {
    return !(lhs == rhs);
  }

  template <typename U>
  std::ostream& operator<<(std::ostream& stream, const Matrix<U>& matrix) {
    const auto [ m, n ] = matrix.getSize();

    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < n; ++j) {
        stream << matrix.getVal(i, j) << " ";
      }

      stream << std::endl;
    }

    return stream;
  }
}

#endif
