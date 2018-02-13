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

  struct Solve {
    enum Method {
      LU,
      Jacobi,
      Thompson
    };
  };

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
          const binaryDual<T>& valMap = binaryDual<T>(
              [](const uint& a, const uint& b) {
                (void) a;
                (void) b;

                return T(0);
              }
          )
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
      template <typename U>
      Matrix(const Matrix<U>& another);

      /**
       * Swap two rows
       * @param fst First index to swap
       * @param snd Second index to swap
       */
      void swapRows(const uint& fst, const uint& snd);

      /**
       * Swap two cols
       * @param fst First index to swap
       * @param snd Second index to swap
       */
      void swapCols(const uint& fst, const uint& snd);

      /**
       * Multiply a row by a scalar
       * @param idx Index of row
       * @param scalar Scalar to multiply
       */
      template <typename U>
      void multiplyRow(const uint& idx, const U& scalar);

      /**
       * Multiply a col by a scalar
       * @param idx Index of col
       * @param scalar Scalar to multiply
       */
      template <typename U>
      void multiplyCol(const uint& idx, const U& scalar);

      /**
       * Add one row to another and scale the first one before adding it
       * @param fst Index of first row
       * @param snd Index of second row
       * @param scalar Multiplier value (by first row)
       */
      template <typename U>
      void addRow(const uint& fst, const uint& snd, const U& scalar = 1);

      /**
       * Add one row to another and scale the first one before adding it
       * @param fst Index of first col
       * @param snd Index of second col
       * @param scalar Multiplier value (by first col)
       */
      template <typename U>
      void addCol(const uint& fst, const uint& snd, const U& scalar = 1);

      /**
       * Get size
       * @return Tuple of the size
       */
      std::tuple<uint, uint> getSize() const;

      /**
       * Determine whether or not matrix is square
       * @return True if square
       */
      bool isSquare() const;

      /**
       * Determine whether or not matrix is diagonal
       * @return True if diagonal
       */
      bool isDiagonal() const;

      /**
       * Determine whether or not matrix is n-diagonal
       * @return True if diagonal
       */
      bool isNDiagonal(const uint& n) const;

      /**
       * Determine whether or not matrix is diagonally dominant
       * @return True if diagonally dominant
       */
      bool isDiagDom() const;

      /**
       * Get the value at the ith row and jth column
       * @param i Row number
       * @param j Column number
       * @return Value stored at Matrix[i][j]
       */
      T getVal(const uint& i, const uint& j) const;

      /**
       * Get the diagonal as a vector
       * @param i Row number
       * @param j Column number
       * @return Value stored at Matrix[i][j]
       */
      std::vector<T> getDiag() const;

      /**
       * Return a matrix of lower triangular entries
       * @param data List of values along lower triangular
       * @return Lower triangular matrix
       */
      Matrix<T> lTriangular() const;

      /**
       * Return a matrix of upper triangular entries
       * @param data List of values along upper triangular
       * @return Upper triangular matrix
       */
      Matrix<T> uTriangular() const;

      /**
       * Generate the LU factorization of this matrix
       * @return A tuple containing the permutation matrix,
       *  lower diagonal, and upper diagonal matrix (respectively)
       */
      std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> LUFactorize() const;

      /**
       * Find the largest eigenvalue for this matrix, and its eigenvector
       * @param nIter Number of iterations to use in power iteration
       * @return A tuple containing the eigenvalue, eigenvector pair
       */
      std::tuple<T, Matrix<T>> largestEigenpair(const uint& nIter = 100);

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

      /* ===== Public Static Methods ===== */

      /**
       * Create a diagonal matrix from a vector
       * @param data List of values along diagonal
       * @return Diagonal matrix
       */
      static Matrix<T> diagonal(const std::vector<T>& data);

      /**
       * Create an identity matrix
       * @param m Size of identity
       * @return Identity matrix
       */
      static Matrix<T> identity(const uint& m);

      /**
       * Create a Hilbert Matrix
       * @param m Size of matrix
       * @return Hilbert matrix of size m
       */
      static Matrix<T> hilbert(const uint& m);

      /**
       * Solve a linear system
       * @param A Linear operator being applied in Ax = b
       * @param b Column vector equal to result
       * @param method Method to solve system with
       * @return Column vector with solution to system
       */
      static Matrix<T> solve(
          const Matrix<T>& A,
          const Matrix<T>& b,
          const Solve::Method& method = Solve::LU
      );

      /**
       * Find the norm of a vector (row or column)
       * @param v Column or row matrix
       * @param n The order of the norm
       */
      static T vNorm(const Matrix<T>& v, const uint& n = 2);

      /**
       * Find the norm of a matrix
       * @param A Matrix to compute the norm of
       * @param n The order of the norm
       */
      static T mNorm(const Matrix<T>& A, const uint& n = 2);

      /**
       * Generate a vector of finite difference coefficients
       * @param order Order of derivative
       * @param accuracy Order of accuracy
       */
      static std::vector<T> genFDCoeff(
          const uint& order,
          const uint& accuracy
      );

      /**
       * Generate a matrix approximation to a differential
       *  operator with a given order and accuracy over
       *  a mesh
       * @param meshSize Size of mesh being computed over
       * @param order Order of derivative
       * @param accuracy Order of accuracy
       */
      static Matrix<T> genFDMatrix(
          const uint& size,
          const uint& order,
          const uint& accuracy = 2
      );

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
      template <typename U, typename V>
      friend Matrix<U> operator*(const V& lhs, const Matrix<U>& rhs);

      /**
       * Scalar multiply a matrix
       * @param lhs Matrix to scale
       * @param rhs Scalar value
       * @return Their matrix sum
       */
      template <typename U, typename V>
      friend Matrix<U> operator*(const Matrix<U>& lhs, const V& rhs);

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
      /* ===== Private Methods ===== */

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

      /* ===== Private Variables ===== */

      // Size of matrix
      uint _m, _n;

      // Storage for values
      std::vector<std::vector<T>> _matrix;
  };
};

#include "matrix.hpp"
#include "matrixOps.hpp"

#endif
