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

  template <typename U>
  Matrix<U> operator*(const U& lhs, const Matrix<U>& rhs) {
    return rhs.scalarMult(lhs);
  }

  template <typename U>
  Matrix<U> operator*(const Matrix<U>& lhs, const U& rhs) {
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
