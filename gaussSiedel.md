---
math: true
permalink: /gaussSiedel
title: Gauss Siedel
layout: page
---

**Routine Name**: Gauss Siedel Solver

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

Gauss Siedel is another matrix separation iterative technique (like Jacobi) that iterates the following until a stable point is found:

# \\[\vec{x}^{(k+1)} = L^{-1} (\vec{b} - U\vec{x}^{(k)})\\]

Where \\(L, U\\) represent the lower-triangular and strictly upper-triangular elements of \\(A\\) in the system \\(A\vec{x} = \vec{b}\\).

**Input**:

A matrix \\(A\\) and result of a linear system as a column vector \\(\vec{b}\\).

**Output**:

A column vector that represents the solution to the linear system.

**Usage/Example**:

{% highlight C++ %}
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

  Mtx _x({{1},{2},{3},{4}});

  auto b = A * _x;

  auto x = Mtx::solve(A, b, Matrix::Solve::GaussSiedel);

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

{% highlight c++ %}
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
      /** Code Omitted **/
    } else if (method == Solve::GaussSiedel) {
      if (!A.isDiagDom()) {
        throw std::domain_error("Input matrix is not diagonally dominant.");
      }

      Matrix<T> x(m, 1, [](const uint& a, const uint& b) {
        (void) a;
        (void) b;

        return 1;
      });

      for (uint iter = 0; iter < 500; ++iter) {
        // Forward substitution algorithm from Wikipedia
        // https://en.wikipedia.org/wiki/Gauss–Seidel_method#Algorithm
        for (uint i = 0; i < m; ++i) {
          T acc = 0;

          for (uint j = 0; j < m; ++j) {
            if (i != j) {
              acc += A.getVal(i, j) * x.getVal(j, 0);
            }
          }

          const auto mult = 1.0 / (double) A.getVal(i, i);
          const auto bi = b.getVal(i, 0);

          x.setVal(i, 0, mult * (bi - acc));
        }
      }

      return x;
    } else if (method == Solve::Thompson) {
      /** Code Omitted **/
    } else if (method == Solve::LU) {
      /** Code Omitted **/
    }

    return Matrix<T>(m);
  }
{% endhighlight %}
