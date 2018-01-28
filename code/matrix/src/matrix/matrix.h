#ifndef MATRIX_H
#define MATRIX_H

#include <functional>
#include <vector>
#include <tuple>

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
       * Copy constructor
       * @param another Another matrix
       */
      Matrix(const Matrix<T>& another);

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

      /* ===== Operators ===== */

      /**
       * Add two matrices
       * @param lhs First matrix to add
       * @param rhs Second matrix to add
       * @return Their matrix sum
       */
      template <typename U>
      friend Matrix<U> operator+(const Matrix<U>& lhs, const Matrix<U>& rhs);

      /**
       * Subtract two matrices
       * @param lhs First matrix
       * @param rhs Second matrix to subtract from a
       * @return Their matrix sum
       */
      template <typename U>
      friend Matrix<U> operator-(const Matrix<U>& lhs, const Matrix<U>& rhs);

      /**
       * Multiply two matrices
       * @param lhs Scalar value
       * @param rhs Matrix to scale
       * @return Their matrix sum
       */
      template <typename U>
      friend Matrix<U> operator*(const Matrix<U>& lhs, const Matrix<U>& rhs);

      /**
       * Scalar multiply a matrix
       * @param lhs Scalar value
       * @param rhs Matrix to scale
       * @return Their matrix sum
       */
      template <typename U>
      friend Matrix<U> operator*(const U& lhs, const Matrix<U>& rhs);

      /**
       * Scalar multiply a matrix
       * @param lhs Matrix to scale
       * @param rhs Scalar value
       * @return Their matrix sum
       */
      template <typename U>
      friend Matrix<U> operator*(const Matrix<U>& lhs, const U& rhs);

      /**
       * Compare matrix equality
       * @param lhs First matrix to compare
       * @param rhs Second matrix to compare
       * @return True if equal
       */
      template <typename U>
      friend bool operator==(const Matrix<U>& lhs, const Matrix<U>& rhs);

      /**
       * Compare matrix equality
       * @param lhs First matrix to compare
       * @param rhs Second matrix to compare
       * @return True if not equal
       */
      template <typename U>
      friend bool operator!=(const Matrix<U>& lhs, const Matrix<U>& rhs);

      /**
       * Output a matrix to the output stream
       * @param stream The output stream
       * @param matrix The matrix to output
       * @return The stream
       */
      template <typename U>
      friend std::ostream& operator<<(
          std::ostream& stream,
          const Matrix<U>& matrix
      );
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
};

#include "matrix.hpp"
#include "matrixOps.hpp"

#endif
