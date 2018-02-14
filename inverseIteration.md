---
math: true
permalink: /inverseIteration
title: Inverse Power Iteration Method
layout: page
---

**Routine Name**: Inverse Iteration

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

Determine the smallest eigenvalue of a matrix, as well as its condition number.

**Input**:

A square matrix

**Output**:

There are two methods. One returns the smallest eigenvalue of the matrix, the other one returns the condition number.

**Usage/Example**:

By finding the approximate condition number for various sizes of the matrix arising in the second order finite difference method approximation, we can see that the condition number grows larger with the size of the mesh. It also appears that the rate does not slow down either. This may indicate that this matrix is poorly conditioned as the size increases.

{% highlight C++ %}
#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  for (uint i = 0; i < 10; ++i) {
    const auto size = (i + 1) * 10;
    auto A = Mtx::genFDMatrix(size, 2);
    const auto smallEigen = std::get<0>(A.smallestEigenpair());
    const auto condition = A.conditionNumber();

    std::cout << "2nd order FDM matrix smallest eigenvalue for size ";
    std::cout << size << ":" << std::endl;
    std::cout << smallEigen << std::endl;

    std::cout << "2nd order FDM matrix condition number for size ";
    std::cout << size << ":" << std::endl;
    std::cout << condition << std::endl;
  }

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
2nd order FDM matrix smallest eigenvalue for size 10:
0.0810141
2nd order FDM matrix condition number for size 10:
45.4552
2nd order FDM matrix smallest eigenvalue for size 20:
0.0223383
2nd order FDM matrix condition number for size 20:
175.087
2nd order FDM matrix smallest eigenvalue for size 30:
0.0102614
2nd order FDM matrix condition number for size 30:
385.728
2nd order FDM matrix smallest eigenvalue for size 40:
0.0058684
2nd order FDM matrix condition number for size 40:
676.365
2nd order FDM matrix smallest eigenvalue for size 50:
0.00379334
2nd order FDM matrix condition number for size 50:
1046.62
2nd order FDM matrix smallest eigenvalue for size 60:
0.00265182
2nd order FDM matrix condition number for size 60:
1497.17
2nd order FDM matrix smallest eigenvalue for size 70:
0.00195755
2nd order FDM matrix condition number for size 70:
2028.16
2nd order FDM matrix smallest eigenvalue for size 80:
0.00150409
2nd order FDM matrix condition number for size 80:
2639.61
2nd order FDM matrix smallest eigenvalue for size 90:
0.00119172
2nd order FDM matrix condition number for size 90:
3331.51
2nd order FDM matrix smallest eigenvalue for size 100:
0.000967435
2nd order FDM matrix condition number for size 100:
4103.86
{% endhighlight %}

**Implementation/Code:**

Finding the smallest eigenvalue is the same as the iterative method described before to find the largest eigenvalue, only now we consider the action of the inverse of the matrix instead of the matrix itself. The result will be an eigenvector of the inverse, which is also an eigenvector of the matrix proper. We extract the eigenvalue the same way we did with forward iteration, and can now determine the condition number as

# \\[ K := \frac{\lambda_{\text{max}}}{\lambda_{\text{min}}} \\].

{% highlight c++ %}
template <typename T>
std::tuple<T, Matrix<T>> Matrix<T>::smallestEigenpair(const uint& nIter) {
  const auto M = std::get<0>(this->getSize());

  // Starting value
  Matrix<T> b(M, 1, [](const uint& i, const uint& j) {
      (void) i;
      (void) j;
      return 1;
  });

  for (uint i = 0; i < nIter; ++i) {
    const auto x = Matrix<T>::solve(*this, b);
    T mult = 1.0 / Matrix<T>::vNorm(x, 2);
    b = x * mult;
  }

  const auto axe = *this * b;
  const auto eigenVal = Matrix<T>::vNorm(axe, 2) / Matrix<T>::vNorm(b, 2);

  return std::tie(eigenVal, axe);
}

template <typename T>
T Matrix<T>::conditionNumber(const uint& nIter) {
  const auto bigEigen = std::get<0>(this->largestEigenpair(nIter));
  const auto smolEigen = std::get<0>(this->smallestEigenpair(nIter));

  return bigEigen / smolEigen;
}
{% endhighlight %}
