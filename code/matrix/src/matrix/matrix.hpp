#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "matrix.h"
#include <iostream>
#include <cmath>
#include <limits>

namespace Matrix {
  /* ===== Utility Functions ===== */
  /**
   * Compute the factorial of a number
   * @param n The number to compute the factorial of
   * @return The factorial
   */
  template <typename T>
  int factorial(uint n) {
    T out = 1;

    for(; n > 1; out *= n--);

    return out;
  }

  template <typename T>
  Matrix<T>::Matrix(
      const uint& m,
      const uint& n,
      const binaryDual<T>& valMap
  ) : _m(m), _n(n) {
    // Initialize 2D matrix to zeroes
    this->_matrix = std::vector<std::vector<T>>(m, std::vector<T>(n, (T) 0));

    // Fill with default funciton
    this->fillWith(valMap);
  }

  template <typename T>
  Matrix<T>::Matrix(const std::vector<std::vector<T>>& init) :
      Matrix(init.size(), init[0].size(), [&](const uint& a, const uint& b) {
        if (init[a].size() == init[0].size()) {
          return init[a][b];
        } else {
          throw std::out_of_range("2D array must be rectangular.");
        }
      }) { };

  template <typename T>
  template <typename U>
  Matrix<T>::Matrix(const Matrix<U>& another) : 
      Matrix(
          std::get<0>(another.getSize()),
          std::get<1>(another.getSize()),
          [&](const uint& a, const uint& b) {
            return (T) another.getVal(a, b);
          }
      ) { }

  /* ===== Public Methods ===== */

  template <typename T>
  void Matrix<T>::swapRows(const uint& fst, const uint& snd) {
    const auto m = std::get<1>(this->getSize());

    if (fst >= 0 && fst < m && snd >= 0 && snd < m) {
      std::swap(this->_matrix[fst], this->_matrix[snd]);
    } else {
      throw std::out_of_range("Indices out of range.");
    }
  }

  template <typename T>
  void Matrix<T>::swapCols(const uint& fst, const uint& snd) {
    const auto [ m, n ] = this->getSize();

    if (fst >= 0 && fst < n && snd >= 0 && snd < n) {
      for (uint i = 0; i < m; ++i) {
        std::swap(this->_matrix[i][fst], this->_matrix[i][snd]);
      }
    } else {
      throw std::out_of_range("Indices out of range.");
    }
  }

  template <typename T>
  template <typename U>
  void Matrix<T>::multiplyRow(const uint& idx, const U& scalar) {
    const auto m = std::get<1>(this->getSize());

    if (idx >= 0 && idx < m) {
      std::transform(
          this->_matrix[idx].begin(),
          this->_matrix[idx].end(),
          this->_matrix[idx].begin(),
          [&](T val) -> T { return val * scalar; }
      );
    } else {
      throw std::out_of_range("Index out of range.");
    }
  }

  template <typename T>
  template <typename U>
  void Matrix<T>::multiplyCol(const uint& idx, const U& scalar) {
    const auto [ m , n ] = this->getSize();

    if (idx >= 0 && idx < n) {
      for (uint i = 0; i < m; ++i) {
        const T val = scalar * this->getVal(i, idx);

        this->setVal(i, idx, val);
      }
    } else {
      throw std::out_of_range("Index out of range.");
    }
  }

  template <typename T>
  template <typename U>
  void Matrix<T>::addRow(const uint& fst, const uint& snd, const U& scalar) {
    const auto [ m , n ] = this->getSize();

    if (fst >= 0 && fst < m && snd >= 0 && snd < m) {
      for (uint i = 0; i < n; ++i) {
        const T val = scalar * this->getVal(fst, i) + this->getVal(snd, i);

        this->setVal(snd, i, val);
      }
    } else {
      throw std::out_of_range("Indices out of range.");
    }
  }

  template <typename T>
  template <typename U>
  void Matrix<T>::addCol(const uint& fst, const uint& snd, const U& scalar) {
    const auto [ m , n ] = this->getSize();

    if (fst >= 0 && fst < n && snd >= 0 && snd < n) {
      for (uint i = 0; i < m; ++i) {
        const T val = scalar * this->getVal(i, fst) + this->getVal(i, snd);

        this->setVal(i, snd, val);
      }
    } else {
      throw std::out_of_range("Indices out of range.");
    }
  }

  template <typename T>
  std::tuple<uint, uint> Matrix<T>::getSize() const {
    return { this->_m, this->_n };
  }

  template <typename T>
  bool Matrix<T>::isSquare() const {
    const auto [ m, n ] = this->getSize();

    return m == n;
  }

  template <typename T>
  bool Matrix<T>::isDiagonal() const {
    return this->isNDiagonal(1);
  }

  template <typename T>
  bool Matrix<T>::isNDiagonal(const uint& n) const {
    if (!this->isSquare()) {
      throw std::domain_error("Matrix cannot be diagonal if square");
    } else if (n % 2 == 0) {
      throw std::domain_error("N-diagonal must have odd n (symmetric).");
    }

    const auto m = std::get<0>(this->getSize());
    const auto base = m - (n - 1) / 2 - 1;

    for (uint col = 0; col < base; ++col) {
      for (uint height = 0; height < base - col; ++height) {
        const auto row = m - height - 1;

        if (this->getVal(row, col) != 0 || this->getVal(col, row) != 0) {
          return false;
        }
      }
    }

    return true;
  }

  template <typename T>
  bool Matrix<T>::isDiagDom() const {
    const auto [ m, n ] = this->getSize();

    for (uint i = 0; i < m; ++i) {
      T sum = 0;

      for (uint j = 0; j < n; ++j) {
        if (i != j) {
          sum += pow(this->getVal(i, j), 2);
        }
      }
      
      if (pow(this->getVal(i, i), 2) < sum) {
        return false;
      }
    }

    return true;
  }

  template <typename T>
  std::vector<T> Matrix<T>::getDiag() const {
    if (this->isSquare()) {
      const auto m = std::get<1>(this->getSize());
      std::vector<T> out;

      for (uint i = 0; i < m; ++i) {
        out.push_back(this->getVal(i, i));
      }

      return out;
    } else {
      throw std::out_of_range("Cannot get diagonal of non-square matrix.");
    }
  }

  template <typename T>
  Matrix<T> Matrix<T>::lTriangular() const {
    const auto [ m, n ] = this->getSize();

    return Matrix<T>(m, n, [&](const uint& a, const uint& b) {
      return a < b ? this->getVal(a, b) : 0;
    });
  }

  template <typename T>
  Matrix<T> Matrix<T>::uTriangular() const {
    const auto [ m, n ] = this->getSize();

    return Matrix<T>(m, n, [&](const uint& a, const uint& b) {
      return a > b ? this->getVal(a, b) : 0;
    });
  }

  template <typename T>
  std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> Matrix<T>::LUFactorize() const {
    if (!this->isSquare()) {
      throw std::out_of_range("Cannot LU factorize non-square matrix.");
    }

    // Initialize outputs
    const auto m = std::get<0>(this->getSize());
    auto P = Matrix<T>::identity(m);
    auto L = Matrix<T>::identity(m);
    auto U = *this;

    // For each column
    for (uint col = 0; col < m; ++col) {
      // Find maximum absolute value in rows below row i
      T max = 0;
      uint swp = col;

      for (uint row = col; row < m; ++row) {
        auto val = U.getVal(row, col);
        auto sqVal = val * val;
        if (sqVal > max) {
          max = sqVal;
          swp = row;
        }
      }

      if (swp != col) {
        // Swap rows with max value to make it a pivot
        U.swapRows(swp, col);
        
        // Store permutation in P matrix
        P.swapRows(swp, col);
      } else if (U.getVal(swp, col) == 0) {
        throw std::domain_error("Matrix is singular. Cannot LU factor.");
      }
      
      // Eliminate all values below pivot
      const auto diag = U.getVal(col, col);

      for (uint row = col + 1; row < m; ++row) {
        const auto mult =  - (U.getVal(row, col) / diag);

        // Eliminate value
        U.addRow(col, row, mult);

        // Add elementary matrix of operation to L
        L.setVal(row, col, -mult);
      }
    }
    
    return std::tuple<Matrix<T>, Matrix<T>, Matrix<T>>(P, L, U);
  }

  template <typename T>
  void Matrix<T>::fillWith(const binaryDual<T>& valMap) {
    const auto [ m , n ] = this->getSize();

    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < n; ++j) {
        this->setVal(i, j, valMap(i, j));
      }
    }
  }

  template <typename T>
  T Matrix<T>::getVal(const uint& i, const uint& j) const {
    if (isInBounds(i, j)) {
      return this->_matrix[i][j];
    } else {
      throw std::out_of_range("Matrix index out of range.");
    }
  }

  template <typename T>
  T Matrix<T>::trace() const {
    const auto m = std::get<1>(this->getSize());
    T sum = 0;

    if (this->isSquare()) {
      for (uint i = 0; i < m; ++i) {
        sum += this->_matrix[i][i];
      }

      return sum;
    } else {
      throw std::domain_error("Matrix must be square to find trace.");
    }
  }

  template <typename T>
  void Matrix<T>::setVal(const uint& i, const uint& j, const T& val) {
    if (isInBounds(i, j)) {
      this->_matrix[i][j] = val;
    } else {
      throw std::out_of_range("Matrix index out of range.");
    }
  }

  template <typename T>
  void Matrix<T>::transpose() {
    const auto [ m, n ] = this->getSize();

    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < i; ++j) {
        std::swap(this->_matrix[i][j], this->_matrix[j][i]);
      }
    }
  }

  /* ===== Public Static Methods ===== */

  template <typename T>
  Matrix<T> Matrix<T>::diagonal(const std::vector<T>& list) {
    return Matrix<T>(
        list.size(),
        list.size(), 
        [&](const uint& a, const uint& b) {
          return (T) (a == b ? list[a] : 0);
        }
    );
  }

  template <typename T>
  Matrix<T> Matrix<T>::identity(const uint& m) {
    return Matrix<T>::diagonal(std::vector<T>(m, (T) 1));
  }

  template <typename T>
  Matrix<T> Matrix<T>::solve(
      const Matrix<T>& A,
      const Matrix<T>& b,
      const Solve::Method& method
  ) {
    // Only work for square matrices
    if (!A.isSquare()) {
      throw std::domain_error("Input matrix must be square.");
    }

    // Get size of matrix
    const auto m = std::get<0>(A.getSize());

    if (method == Solve::Jacobi) {
      if (!A.isDiagDom()) {
        throw std::domain_error("Input matrix is not diagonally dominant.");
      }

      Matrix<T> invD(m, m, [&](const uint& a, const uint& b) -> T {
        if (a == b && A.getVal(a, b) != 0) {
          return 1.0 / A.getVal(a, b);
        } else {
          return 0;
        }
      });

      auto R = A.lTriangular() + A.uTriangular();

      auto x = b;

      for (unsigned int i = 0; i < 500; ++i) {
        x = invD * (b - R * x);
      }

      return x;
    } else if (method == Solve::Thompson) {
      if (!A.isNDiagonal(3)) {
        throw std::domain_error("Thompson method needs tri-diagonal matrix.");
      }

      // Throw diagonals into vectors (store in a vector of diags)
      std::vector<std::vector<T>> diags(3, std::vector<T>());

      diags[0].push_back(0);

      for (uint row = 0; row < m; ++row) {
        for (uint nDiag = 0; nDiag < 3; ++nDiag) {
          const auto pos = row - 1 + nDiag;

          if (pos < m) {
            diags[nDiag].push_back(A.getVal(row, (uint) pos));
          }
        }
      }

      diags[2].push_back(0);

      // Syntactic nicety
      auto _a = diags[0];
      auto _b = diags[1];
      auto _c = diags[2];

      // Copy solution vector
      auto _d = std::vector<T>();
      for (uint i = 0; i < m; ++i) {
        _d.push_back(b.getVal(i, 0));
      }

      // Do elimination via algebra
      _c[0] = _c[0] / _b[0];
      _d[0] = _d[0] / _b[0];
      for (uint i = 1; i < m; ++i) {
        _c[i] = _c[i] / (_b[i] - _a[i] * _c[i - 1]);
        _d[i] = (_d[i] - _a[i] * _d[i - 1]) / (_b[i] - _a[i] * _c[i - 1]);
      }

      // Create the solution
      Matrix<T> x(m, 1);
      x.setVal(m - 1, 0, _d[m - 1]);
      for (int i = m - 2; i >= 0; --i) {
        const uint _i = i;

        x.setVal(_i, 0, _d[_i] - _c[_i] * x.getVal(_i + 1, 0));
      }

      return x;
    } else if (method == Solve::LU) {
      // Factor A into components
      auto [ P, L, U ] = A.LUFactorize();

      // Permute result vector
      P.transpose();
      auto res = P * b;

      /* ===== Solve Ly = res ===== */
      Matrix<T> y(m, 1);

      // For each row (top to bottom)
      for (uint row = 0; row < m; ++row) {
        // Sum up the equation to the left (already solved)
        T leftSum = 0;

        for (int col = row - 1; col >= 0; --col) {
          leftSum += L.getVal(row, col) * y.getVal(col, 0);
        }

        // Solve for the value at that row
        y.setVal(row, 0, (res.getVal(row, 0) - leftSum) / L.getVal(row, row));
      }

      /* ===== Solve Ux = y ===== */
      Matrix<T> x(m, 1);

      // For each row (bottom to top)
      for (int row = m - 1; row >= 0; --row) {
        // Sum up the equation to the right (already solved)
        T rightSum = 0;

        for (uint col = row + 1; col < m; ++col) {
          rightSum += U.getVal(row, col) * x.getVal(col, 0);
        }

        // Solve for the value at that row
        x.setVal(row, 0, (y.getVal(row, 0) - rightSum) / U.getVal(row, row));
      }

      return x;
    }

    return Matrix<T>(m);
  }

  template <typename T>
  T Matrix<T>::vNorm(const Matrix<T>& v, const uint& n) {
    const auto [ M, N ] = v.getSize();
    if (M != 1 && N != 1 && M != N) {
      throw std::domain_error("Need a row or column vector for vector norm.");
    }

    const auto isRow = (M == 1);
    const auto size = isRow ? N : M;

    T sum = 0;
    T max = -std::numeric_limits<T>::max();

    for (uint i = 0; i < size; ++i) {
      const auto val = isRow ? v.getVal(0, i) : v.getVal(i, 0);

      max = (val > max) ? val : max;
      sum += pow(std::abs(val), n);
    }

    if (n == std::numeric_limits<uint>::max()) {
      // Infinity norm
      return max;
    } else {
      return pow(sum, 1.0 / n);
    }
  }

  template <typename T>
  std::vector<T> Matrix<T>::genFDCoeff(const uint& order, const uint& accuracy) {
    // Determine the amount of coefficients needed
    auto size = 2 * std::floor((order + 1) / 2) - 1 + accuracy;

    // Determine max absolute index count
    const auto p = (size - 1) / 2;

    // Initialize the taylor system matrix
    const auto A = Matrix<T>(size, size, [&](const uint& a, const uint& b) {
        // Return taylor coefficient at this value
        return std::pow((-p + b), a);
    });

    // Initialize result vector
    Matrix<T> b(size, 1);
    b.setVal(order, 0, factorial<uint>(order));

    const auto x = Matrix<T>::solve(A, b);
    std::vector<T> coeffs;
    for (uint i = 0; i < size; ++i) {
      coeffs.push_back(x.getVal(i, 0));
    }

    return coeffs;
  }

  template <typename T>
  Matrix<T> Matrix<T>::genFDMatrix(
      const uint& size,
      const uint& order,
      const uint& accuracy
  ) {
    const auto coeffs = Matrix<T>::genFDCoeff(order, accuracy);

    return Matrix<T>(size, size, [&](const uint& row, const uint& col) -> T {
        const int start = row - 1;
        const int end = start + coeffs.size() - 1;

        if ((int) col >= start && (int) col <= end) {
          return coeffs[(int) col - row + 1];
        } else {
          return 0;
        }
    });
  }

  /* ===== Private Methods ===== */

  template <typename T>
  Matrix<T> Matrix<T>::add(const Matrix<T>& another) const {
    // Determine compatibility
    const auto [ m, n ] = this->getSize();
    const auto [ M, N ] = another.getSize();

    if (m == M && n == N) {
      return Matrix<T>(m, n, [&](const uint& a, const uint& b) {
        return this->getVal(a, b) + another.getVal(a, b);
      });
    } else {
      throw std::out_of_range("Can not be added, wrong dimensions.");
    }
  }

  template <typename T>
  Matrix<T> Matrix<T>::subtract(const Matrix<T>& another) const {
    // Determine compatibility
    const auto [ m, n ] = this->getSize();
    const auto [ M, N ] = another.getSize();

    if (m == M && n == N) {
      return Matrix<T>(m, n, [&](const uint& a, const uint& b) {
        return this->getVal(a, b) - another.getVal(a, b);
      });
    } else {
      throw std::out_of_range("Can not be subtracted, wrong dimensions.");
    }
  }

  template <typename T>
  Matrix<T> Matrix<T>::scalarMult(const T& scalar) const {
    const auto [ m, n ] = this->getSize();

    return Matrix<T>(m, n, [&](const uint& a, const uint& b) {
      return scalar * this->getVal(a, b);
    });
  }

  template <typename T>
  Matrix<T> Matrix<T>::multiply(const Matrix<T>& another) const {
    // Determine compatibility
    const auto [ m, n ] = this->getSize();
    const auto [ M, N ] = another.getSize();

    if (n == M) {
      auto out = Matrix<T>(m, N);

      for (uint i = 0; i < m; ++i) {
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
    const auto [ m, n ] = this->getSize();
    const auto [ M, N ] = another.getSize();

    if (m == M && n == n) {
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
    const auto [ m, n ] = this->getSize();

    return (i >= 0) && (i < m) && (j >= 0) && (j < n);
  }
};

#endif
