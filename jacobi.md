---
math: true
permalink: /jacobi
title: Jacobi Iteration
layout: page
---

**Routine Name**: Jacobi Iteration

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This algorithm will solve a linear system such that the linear operator is diagonally dominant.

**Input**:

A diagonally dominant matrix \\( A \\) and a column vector \\( b \\) compatible with it.

**Output**:

The solution vector \\( x \\) such that \\( Ax = b \\).

**Usage/Example**:

Here I pre-generate a dummy vector, then operate on it with a tri-diagonal matrix \\( A \\) to get my \\( b \\) vector. This means that the solution should be very close to the dummy vector I supplied. As you can see from the output, it works well.

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

  auto x = Mtx::solve(A, b, Matrix::Solve::Jacobi);

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

Jacobi iteration is the continued iteration of:

# \\[ x_{k+1} = D^{-1} (b - Rx_k) \\]

Where \\( D \\) is the diagonal of \\( A \\), \\( R \\) is \\( A \\) without the diagonal. Luckily, the inverse of a diagonal matrix can be found in \\( O(n) \\) by inverting every element in the diagonal (or just do nothing if the element is \\( 0 \\) ). Since matrix multiplication of a column vector takes \\( O(n^2) \\) operations, this algorithm runs in

# \\[ O(n^2) \\]

{% highlight c++ %}
template <typename T>
Matrix<T> Matrix<T>::solve(
    const Matrix<T>& A,
    const Matrix<T>& b
) {
  // Only work for square matrices
  if (!A.isSquare()) {
    throw std::domain_error("Input matrix must be square.");
  }

  // Get size of matrix
  const auto m = std::get<0>(A.getSize());

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
}
{% endhighlight %}
