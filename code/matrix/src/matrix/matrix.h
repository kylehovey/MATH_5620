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
  const binaryDual<T> zero = binaryDual<T>(
      [](const uint& a, const uint& b) {
        (void) a;
        (void) b;

        return (T) 0;
      }
  );

  /* ===== Class Definition ===== */

  template <typename T>
  class Matrix {
    public:
      /**
       * Construct using size and fill function
       * @param m Number of rows
       * @param n Number of columns
       * @param valMap Function to fill matrix with initial values
       */
      Matrix(
          const uint& m = 3,
          const uint& n = 3,
          const binaryDual<T>& valMap = zero<T>
      );

      /**
       * Construct using a 2D matrix
       * @param init 2D vector to build matrix from
       */
      Matrix(const std::vector<std::vector<T>>& init);

      /**
       * Get size
       * @return Tuple of the size
       */
      std::tuple<uint, uint> getSize() const;

      /**
       * Get the value at the ith row and jth column
       * @param i Row number
       * @param j Column number
       * @return Value stored at Matrix[i][j]
       */
      T getVal(const uint& i, const uint& j) const;

      /**
       * Get the trace of this matrix
       * @return The trace
       */
      T trace() const;

      /**
       * Set the value at the ith row and jth column
       * @param i Row number
       * @param j Column number
       * @param val Value to set
       */
      void setVal(const uint& i, const uint& j, const T& val);

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
      /**
       * Determine whether or not values are within bounds
       * @param i Row number
       * @param j Column number
       * @return True if in range, false if not
       */
      bool isInBounds(const uint& i, const uint& j) const;

      /**
       * Add another matrix to this one and return the result
       * @param another Another matrix
       * @return This matrix + another
       */
      Matrix<T> add(const Matrix<T>& another) const;

      /**
       * Subtract another matrix from this one and return the result
       * @param another Another matrix
       * @return This matrix - another
       */
      Matrix<T> subtract(const Matrix<T>& another) const;

      /**
       * Multiply this matrix by a scalar and return the result
       * @param scalar Value to multiply this matrix by
       * @return Scalar * this matrix
       */
      Matrix<T> scalarMult(const T& scalar) const;

      /**
       * Multiply another matrix by this one and return the result
       * @param another Another matrix
       * @return This matrix * another
       */
      Matrix<T> multiply(const Matrix<T>& another) const;

      /**
       * Determine whether or not another matrix is equal to this one
       * @param another Another matrix
       * @return True if matrix is equal
       */
      bool isEqualTo(const Matrix<T>& another) const;

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

  template <typename T>
  Matrix<T>::Matrix(const std::vector<std::vector<T>>& init) :
      Matrix(init.size(), init[0].size(), [&](const uint& a, const uint& b) {
            return init[a][b];
      }) { };

  /* ===== Public Methods ===== */

  template <typename T>
  std::tuple<uint, uint> Matrix<T>::getSize() const {
    return { this->m, this->n };
  }

  template <typename T>
  void Matrix<T>::fillWith(const binaryDual<T>& valMap) {
    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < n; ++j) {
        this->setVal(i, j, valMap(i, j));
      }
    }
  }

  template <typename T>
  T Matrix<T>::getVal(const uint& i, const uint& j) const {
    if (isInBounds(i, j)) {
      return this->matrix[i][j];
    } else {
      throw std::out_of_range("Matrix index out of range.");
    }
  }

  template <typename T>
  T Matrix<T>::trace() const {
    T sum = 0;

    if (this->m == this->n) {
      for (uint i = 0; i < this->m; ++i) {
        sum += this->matrix[i][i];
      }

      return sum;
    } else {
      throw std::domain_error("Matrix must be square to find trace.");
    }
  }

  template <typename T>
  void Matrix<T>::setVal(const uint& i, const uint& j, const T& val) {
    if (isInBounds(i, j)) {
      this->matrix[i][j] = val;
    } else {
      throw std::out_of_range("Matrix index out of range.");
    }
  }

  template <typename T>
  void Matrix<T>::transpose() {
    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < i; ++j) {
        std::swap(this->matrix[i][j], this->matrix[j][i]);
      }
    }
  }

  /* ===== Private Methods ===== */

  template <typename T>
  Matrix<T> Matrix<T>::add(const Matrix<T>& another) const {
    // Determine compatibility
    const auto [ M, N ] = another.getSize();

    if (this->m == M && this->m == M) {
      return Matrix<T>(this->m, this->n, [&](const uint& a, const uint& b) {
        return this->getVal(a, b) + another.getVal(a, b);
      });
    } else {
      throw std::out_of_range("Can not be added, wrong dimensions.");
    }
  }

  template <typename T>
  Matrix<T> Matrix<T>::subtract(const Matrix<T>& another) const {
    // Determine compatibility
    const auto [ M, N ] = another.getSize();

    if (this->m == M && this->m == M) {
      return Matrix<T>(this->m, this->n, [&](const uint& a, const uint& b) {
        return this->getVal(a, b) - another.getVal(a, b);
      });
    } else {
      throw std::out_of_range("Can not be subtracted, wrong dimensions.");
    }
  }

  template <typename T>
  Matrix<T> Matrix<T>::scalarMult(const T& scalar) const {
    return Matrix<T>(this->m, this->n, [&](const uint& a, const uint& b) {
      return scalar * this->getVal(a, b);
    });
  }

  template <typename T>
  Matrix<T> Matrix<T>::multiply(const Matrix<T>& another) const {
    // Determine compatibility
    const auto [ M, N ] = another.getSize();

    if (this->n == M) {
      auto out = Matrix<T>(this->m, N);

      for (uint i = 0; i < this->m; ++i) {
        for (uint j = 0; j < N; ++j) {
          T sum = 0;

          for (uint k = 0; k < M; ++k) {
            sum += this->getVal(i, k) * another.getVal(k, j);
          }

          out.setVal(i, j, sum);
        }
      }

      return out;
    } else {
      throw std::out_of_range("Matrices can not be multiplied.");
    }
  }

  template <typename T>
  bool Matrix<T>::isEqualTo(const Matrix<T>& another) const {
    const auto [ M, N ] = another.getSize();

    if (this->m == M && this->n == n) {
      for (uint i = 0; i < M; i++) {
        for (uint j = 0; j < N; j++) {
          if (this->getVal(i, j) != another.getVal(i, j)) {
            return false;
          }
        }
      }
    } else {
      return false;
    }

    return true;
  }

  template <typename T>
  bool Matrix<T>::isInBounds(const uint& i, const uint& j) const {
    return (i >= 0) && (i < this->m) && (j >= 0) && (j < this->n);
  }
};

#endif
