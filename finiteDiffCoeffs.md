---
math: true
permalink: /finiteDiffCoeffs
title: Hardware Epsilon
layout: page
---

**Routine Name**: genFDCoeff

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This method generates the coefficients required for a finite difference method approximating a derivative of order `order` and of accuracy `accuracy`.

**Input**:

The order and accuracy desired for the FDM coefficients.

**Output**:

A vector containing the finite difference coefficients.

**Usage/Example**:

{% highlight C++ %}
#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  std::cout << "1st derivative (m = 1) with order 4 accuracy (n = 4)\n";
  const auto coeffs = Mtx::genFDCoeff(2, 2);

  for (auto i = 0u; i < coeffs.size(); ++i) {
    std::cout << coeffs[i] << std::endl;
  }

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
1st derivative (m = 1) with order 2 accuracy (n = 2)
1
-2
1
{% endhighlight %}

**Implementation/Code:**

{% highlight c++ %}
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
{% endhighlight %}
