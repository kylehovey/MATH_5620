#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "matrix.h"

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
        std::swap(this->matrix[i][j], this->matrix[j][i]);
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
