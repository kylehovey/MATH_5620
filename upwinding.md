---
math: true
permalink: /upwinding
title: Upwinding Method
layout: page
---

**Routine Name**: Upwinding Method

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This method is used to solve partial differential equations that result in advection, such as

# \\[ \frac{\partial u}{\partial t} = a \frac{\partial u}{\partial x} \\]

Equations such as these are very easy to solve for (the initial condition just propagates forward/backward in time), but it is clear that finding a finite difference method to model the equation is less than trivial.

**Input**:

* Domain in space
* Domain in time
* Time step
* Space step
* Constant \\(a\\) found in PDE
* Boundary conditions \\( \eta(x) \\) at \\(t_0\\)

**Output**:

An image where space is represented vertically and time horizontally from left to right.

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include "./advection/advection.h"

int main() {
  const auto spaceDomain = std::make_tuple(0.0, 1.0);
  const auto timeDomain = std::make_tuple(0.0, 1.0);
  const auto dt = 1e-3;
  const auto dx = 1e-3;
  const auto [ A, B ] = spaceDomain;
  const auto size = std::round((B - A) / dx);
  const auto a = 0.7;
  const auto eta = [](const double& x) -> double {
    return (x >= 0.3 && x <= 0.6) ? 100 : 0;
  };

  Advection::testUpwinding<double>(
      size,
      spaceDomain,
      timeDomain,
      dx,
      dt,
      eta,
      a
  );

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

Notice the diffusion along the edges of the wave-front. This is due to the jump discontinuity found in the initial conditions and can be improved by adding a second-order diffusive property (as seen in the next Lax-Wendroff example).

![Upwinding Method](/MATH_5620/images/upWinding.png)

**Implementation/Code:**

{% highlight C++ %}
/**
 * Test the upwinding method for the PDE:
 *  u_t + a u_x = 0; boundary conditions: 0
 * @param size Size of mesh
 * @param spaceDomain Interval in space to evaluate over
 * @param timeDomain Interval in time to evaluate over
 * @param spaceStep Space step quantity
 * @param timeStep Time step quantity
 * @param eta Initial conditions along t = t_0
 * @param a Value seen in advection equation
 */
template<typename T>
void testUpwinding(
    const uint size,
    const interval<T> spaceDomain,
    const interval<T> timeDomain,
    const T& spaceStep,
    const T& timeStep,
    const monad<T>& eta,
    const T& a
) {
  // Solution vector
  Matrix::Matrix<T> U(size - 2, 1, [=](const uint& a, const uint& b) -> T {
      (void) b;

      const auto A = std::get<0>(spaceDomain);
      const auto x = A + (a + 1) * spaceStep;

      return eta(x);
  });

  // Time-step operator
  const auto coeffs = (std::vector<T>) {
    a * spaceStep / timeStep,
    1 - a * spaceStep / timeStep
  };
  const auto A = Matrix::Matrix<T>::genNDiag(
      size - 2, 
      coeffs, 
      -1
  );

  // Evolve system in time
  std::vector<Matrix::Matrix<T>> solns = { U };
  const auto [ t_0, t_f ] = timeDomain;
  for (auto t = t_0; t < t_f; t += timeStep) {
    U = A * U;
    solns.push_back(U);
  }

  const Matrix::Matrix<T> wholePicture(
      size - 2,
      solns.size(),
      [=](const uint& a, const uint& b) -> double {
        return solns[b].getVal(a, 0);
      }
  );

  // Output
  Image::ImageWriter::matrixHeatmap("./upWinding.ppm", wholePicture, 1000);
}
{% endhighlight %}
