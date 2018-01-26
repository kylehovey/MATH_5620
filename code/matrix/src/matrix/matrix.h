#ifndef MATRIX_H
#define MATRIX_H

#include <functional>
#include <vector>
#include <iostream>

namespace Matrix {
  /* ===== Namespace Typedef and Variable Definition ===== */
  using uint = unsigned int;

  template <typename T>
  using binaryDual = std::function<T(const uint&, const uint&)>;

  template <typename T>
  const binaryDual<T> kronecker = binaryDual<T>(
      [](const uint& a, const uint& b) {
        return (T) ((a == b) ? 1 : 0);
      }
  );

  /* ===== Class Definition ===== */

  template <typename T>
  class Matrix {
    public:
      /**
       * @param m Number of rows
       * @param n Number of columns
       * @param valMap Function to fill matrix with initial values
       */
      Matrix(
          const uint& m = 3,
          const uint& n = 3,
          const binaryDual<T>& valMap = kronecker<T>
      );

      /**
       * Get the value at the ith row and jth column
       * @param i Row number
       * @param j Column number
       * @return Value stored at Matrix[i][j]
       */
      T getVal(const uint& i, const uint& j) const;

      /**
       * Fill the matrix with values determined by indices
       * @param valMap Function to fill matrix with values
       */
      void fillWith(const binaryDual<T>& valMap);

      /**
       * Transpose this matrix
       */
      void transpose();
    private:
      // Size of matrix
      uint m, n;

      // Storage for values
      std::vector<std::vector<T>> matrix;
  };

  /* ===== Constructor Definition ===== */

  template <typename T>
  Matrix<T>::Matrix(
      const uint& m,
      const uint& n,
      const binaryDual<T>& valMap
  ) : m(m), n(n) {
    // Initialize 2D matrix to zeroes
    this->matrix = std::vector<std::vector<T>>(m, std::vector<T>(n, (T) 0));

    // Fill with default funciton
    this->fillWith(valMap);
  }

  /* ===== Public Methods ===== */

  template <typename T>
  void Matrix<T>::fillWith(const binaryDual<T>& valMap) {
    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < n; ++j) {
        this->matrix[i][j] = valMap(i, j);
      }
    }
  }

  template <typename T>
  T Matrix<T>::getVal(const uint& i, const uint& j) const {
    return this->matrix[i][j];
  }

  template <typename T>
  void Matrix<T>::transpose() {
    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < i; ++j) {
        std::swap(this->matrix[i][j], this->matrix[j][i]);
      }
    }
  }
};

#endif
