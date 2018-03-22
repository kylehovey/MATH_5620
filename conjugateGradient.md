---
math: true
permalink: /conjugateGradient
title: Conjugate Gradient
layout: page
---

**Routine Name**: Conjugate Gradient Method for Solving Linear Systems

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

The Conjugate Gradient method utilizes the vector gradient of a generalized potential function for a linear system \\(A\vec{x} = \vec{b}\\). This method works well for elliptic systems as we can guarantee a global minima. We reach successive minimizers along the direction of the instantaneous gradient of our potential function by finding the length we need to step in that direction to reach a valley.

This method resulted in only a few (usually less than 10) iterations to completion, whereas matrix separation techniques I have used before took hundreds of iterations on average.

**Input**:

The matrix and the output of the linear system.

**Output**:

A column vector as the solution of the linear system

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

  Mtx _x({ {1},{2},{3},{4} });

  auto b = A * _x;

  auto x = Mtx::solve(A, b, Matrix::Solve::ConjugateGradient);

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
    /** Code Omitted **/
  } else if (method == Solve::ConjugateGradient) {
    auto x = b;
    auto r = b - (A * x);
    auto p = r;

    for (uint i = 0; i < 500; ++i) {
      const double alpha =
        Matrix::innerProduct(r, r) / Matrix::innerProduct(p, A * p);

      x = x + (alpha * p);
      const auto lastR = r;
      r = r - (alpha * A * p);

      if (Matrix::vNorm(r - lastR, 2) < 0.001) { break; }

      const double beta =
        Matrix::innerProduct(r, r) / Matrix::innerProduct(lastR, lastR);

      p = r + (beta * p);
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
