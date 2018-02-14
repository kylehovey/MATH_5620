---
math: true
permalink: /powerIteration
title: Power Iteration Method
layout: page
---

**Routine Name**: Power Iteration

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This method finds the eigenvalue with largest absolute magnitude assuming a strict ordering.

**Input**:

A square matrix.

**Output**:

A tuple containing the largest eigenvalue and its corresponding eigenvector.

**Usage/Example**:

{% highlight C++ %}
#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  auto A = Mtx::hilbert(5);

  const auto [ bigEigen, x ] = A.largestEigenpair();

  std::cout << "A\n";
  std::cout << A << std::endl;

  std::cout << "Largest Eigenvalue\n";
  std::cout << bigEigen << std::endl;
  std::cout << std::endl;
  std::cout << "x vector\n";
  std::cout << x << std::endl;
  std::cout << "A * x\n";
  std::cout << A * x << std::endl;
  std::cout << "Largest Eigenvalue * x\n";
  std::cout << bigEigen * x << std::endl;

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
A
1 0.5 0.333333 0.25 0.2
0.5 0.333333 0.25 0.2 0.166667
0.333333 0.25 0.2 0.166667 0.142857
0.25 0.2 0.166667 0.142857 0.125
0.2 0.166667 0.142857 0.125 0.111111

Largest Eigenvalue
1.56705

x vector
1.20327
0.698577
0.503929
0.397152
0.328803

A * x
1.88558
1.09471
0.789683
0.622357
0.515251

Largest Eigenvalue * x
1.88558
1.09471
0.789683
0.622357
0.515251
{% endhighlight %}

**Implementation/Code:**

This uses the forward iterative method

# \\[ \vec{v}_{k+1} = \frac{A \vec{v}_k}{\left\lVert A \vec{v}_k \right\rVert_2} \\]

{% highlight c++ %}
template <typename T>
std::tuple<T, Matrix<T>> Matrix<T>::largestEigenpair(const uint& nIter) {
  const auto M = std::get<0>(this->getSize());

  // Starting value
  Matrix<T> b(M, 1, [](const uint& i, const uint& j) {
      (void) i;
      (void) j;
      return 1;
  });

  for (uint i = 0; i < nIter; ++i) {
    const auto x = *this * b;
    T mult = 1.0 / Matrix<T>::vNorm(x, 2);
    b = mult * x;
  }

  const auto axe = *this * b;
  const auto eigenVal = Matrix<T>::vNorm(axe, 2) / Matrix<T>::vNorm(b, 2);

  return std::tie(eigenVal, axe);
}
{% endhighlight %}
