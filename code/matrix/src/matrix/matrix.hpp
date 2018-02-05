#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "matrix.h"
#include <iostream>

namespace Matrix {
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
    const auto [ m, n ] = this->getSize();

    for (uint i = 0; i < m; ++i) {
      for (uint j = 0; j < n; ++j) {
        if (i != j && this->getVal(i, j) != 0) {
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
          sum += this->getVal(i, j);
        }
      }
      
      if (this->getVal(i, i) < sum) {
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
    // Get size of matrix
    const auto m = std::get<0>(A.getSize());

    if (method == Solve::Jacobi) {
    } else if (method == Solve::Thompson) {
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

    // TODO Dummy return
    return Matrix<T>(m);
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
