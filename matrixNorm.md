---
math: true
permalink: /matrixNorm
title: Matrix Norms
layout: page
---

**Routine Name**: Matrix Norm

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This method computes the matrix \\(1\\)-norm and \\(\infty\\)-norm. The former is equivalent to the maximum column sum, and the latter is equivalent to the maximum row sum.

**Input**:

A square matrix.

**Output**:

The \\(1\\)-norm or \\(\infty\\)-norm of that matrix.

**Usage/Example**:

{% highlight C++ %}
#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <limits>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx A({
      {1, 5, 1, 0},
      {2, 0, 0, 1},
      {3, 8, 5, 2},
      {4, 9, 2, 1}
  });

  auto oneNorm = Mtx::mNorm(A, 1);
  auto infNorm = Mtx::mNorm(A, std::numeric_limits<uint>::max());

  std::cout << "1-norm of matrix:" << std::endl;
  std::cout << oneNorm << std::endl;
  std::cout << "Infinity-norm of matrix:" << std::endl;
  std::cout << infNorm << std::endl;

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
1-norm of matrix:
22
Infinity-norm of matrix:
18
{% endhighlight %}

**Implementation/Code:**

{% highlight c++ %}
template <typename T>
T Matrix<T>::mNorm(const Matrix<T>& A, const uint& n) {
  const auto [ M, N ] = A.getSize();
  T max = -std::numeric_limits<uint>::max();

  if (n == std::numeric_limits<uint>::max()) {
    // Max row-sum
    for (uint row = 0; row < M; ++row) {
      T sum = 0;

      for (uint col = 0; col < N; ++col) {
        sum += A.getVal(row, col);
      }

      max = sum > max ? sum : max;
    }

    return max;
  } else if (n == 1) {
    // Max col-sum
    for (uint col = 0; col < N; ++col) {
      T sum = 0;

      for (uint row = 0; row < N; ++row) {
        sum += A.getVal(row, col);
      }

      max = sum > max ? sum : max;
    }

    return max;
  } else {
    throw std::domain_error("Matrix norm order not implemented.");
  }
}
{% endhighlight %}
