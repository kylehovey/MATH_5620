---
math: true
permalink: /lu
title: LU Factorization
layout: page
---

**Routine Name**: LU Factorization

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This algorithm will solve a linear system such that the linear operator is non-singular and square.

**Input**:

A square matrix \\( A \\) and a column vector \\( b \\) compatible with it.

**Output**:

The solution vector \\( x \\) such that \\( Ax = b \\).

**Usage/Example**:

Here I pre-generate a dummy vector, then operate on it with a tri-diagonal matrix \\( A \\) to get my \\( b \\) vector. This means that the solution should be very close to the dummy vector I supplied. As you can see from the output, it works well.

{% highlight c++ %}
#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx A({
      { -2, 1, 0, 0 },
      { 1, -2, 1, 0 },
      { 0, 1, -2, 1 },
      { 0, 0, 1, -2 }
  });

  Mtx _x({ {1},{2},{3},{4} });

  auto b = A * _x;

  auto x = Mtx::solve(A, b, Matrix::Solve::LU);

  std::cout << x << std::endl;

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
1
2
3
4
{% endhighlight %}

**Implementation/Code:**

LU factorization requires that you do Gaussian Elimination, which requires \\( O(n^3) \\) operations. The only difference here is that I keep track of the operations required and assemble the \\( L \\) and \\( U \\) matrix (as well as the permutation matrix \\( P \\) ) so that they may be returned at the end. Because of this equivalence, this algorithm runs in:

# \\[ O(n^3) \\]

{% highlight c++ %}
template <typename T>
Matrix<T> Matrix<T>::solve(
    const Matrix<T>& A,
    const Matrix<T>& b
) {
  // Factor A into components
  auto [ P, L, U ] = A.LUFactorize();

  // Permute result vector
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
{% endhighlight %}

Where LU factorization is done by:

{% highlight c++ %}
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
{% endhighlight %}
