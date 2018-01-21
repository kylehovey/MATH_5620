---
math: true
permalink: /springmass
title: Spring Mass Differential Equation
layout: page
---

**Routine Name**: Spring Mass Equation

**Author**: Kyle Hovey

**Language**: C++
**Description/Purpose**: 
A spring-mass system with mass \\( m \\), damping coefficient \\( \gamma \\), and spring constant \\( k \\) is described by the differential equation

# \\[ m y\'\'(t) + \gamma y\'(t) + k y(t) = f(t) \\]

where \\( f(t) \\) is defined to be the driving function for the system. This equation is derived from a sum of forces, so \\( f(t) \\) is equal to the force left over when all other forces cancel out.

**Input**:
* Mass `m`
* Damping Coefficient `gamma`
* Spring Constant `k`
* Initial Conditions \\( y_o, y'_o \\)

**Output**:

A function that takes time as input and returns the position of the mass at that time.

**Usage/Example**:

To use this function, you have to include it in the file you need to call it from.

{% highlight C++ %}
#include "springmass/springmass.h"
#include <limits>
#include <iostream>

int main() {
  // Constants and initial values
  const double
      m = 10,
      gamma = 1,
      k = 1,
      yo = 0,
      dyo = 1;

  // Generate the analytic solution
  auto y = genSpringMass<double>(yo, dyo, m, gamma, k);

  // Output position at time t for 3 seconds (damped oscillation)
  for (double t = 0; t < 3; t += 0.1) {
    std::cout << y(t) << std::endl;
  }

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
0
0.493385
0.967246
1.4126
1.82122
2.18574
2.49983
2.75829
2.9571
3.09352
3.16608
3.17462
3.12024
3.00524
2.83312
2.6084
2.33656
2.02392
1.67746
1.30469
0.913521
0.512036
0.108378
-0.289429
-0.673638
-1.03692
-1.37251
-1.67433
-1.93707
-2.15631
{% endhighlight %}

**Implementation/Code:**

The spring-mass differential equation is a second degree linear equation, which will have two roots to its corresponding characteristic equation. These roots should be allowed to have complex domain as a spring mass system is one that often tends to be oscillatory. The complimentary solution will have the form:

# \\[ y_c(t) = C_1 e^{\lambda_1 t} + C_2 e^{\lambda_2 t} \\]

The particular solution would involve quite a bit of integration, and most likely some clever use of the Wronskian determinant via the method of Variation of Parameters. As we are not yet equipped to compute such integration in code in this class, I have left the particular solution as a stub to be filled out later.

Solving for the analytic particular solution to the homogenous part of the spring-mass differential equation, we can eliminate \\( C_1 \\) and \\( C_2 \\) by algebraicly determining their equivalents with respect to \\( y_o \\) and \\( y'_o \\).

All this leaves for programming is determining the eigenvalues and leading coefficients, then to construct a function to determine the position of the mass at time `t`:

{% highlight c++ %}
#include <functional>
#include <complex>
#include <iostream>

// We will be using this type a lot, make an alias for it
template <typename T>
using endomorphism = std::function<T(T)>;

// Make an alias for the zero function
template <typename T>
const endomorphism<T> zeroFunction = endomorphism<T>(
  /**
   * @param t Temporal input
   * @return Zero cast into the template type
   */
  [](T t) {
    (void) t;
    return (T) 0;
  }
);

template <typename T>
/**
 * Generate an analytic solution to a spring mass system that gives
 * the position of the mass at any time t.
 * @param yo Initial position (t = 0)
 * @param dyo Initial derivative (t = 0)
 * @param m The mass
 * @param gamma Damping constant (air resistance, etc...)
 * @param k Spring constant
 * @param endomorphism f Driving function (non-homogenous part)
 * @return Solution to spring-mass equation
 */
endomorphism<T> genSpringMass(
    T yo = 0,
    T dyo = 0,
    T m = 1,
    T gamma = 1,
    T k = 1,
    endomorphism<T> f = zeroFunction<T>
) {
  // Compute the eigenvalues and allow for complex solutions
  auto desc = std::sqrt(
      (std::complex<T>) (std::pow(gamma, 2) - 4 * m * k)
  );

  std::complex<T> eigens[2] = {
    (-gamma + desc) / std::complex<double>(2 * m),
    (-gamma - desc) / std::complex<double>(2 * m)
  };

  // Determine the leading coefficients of the solution
  auto denom = eigens[1] - eigens[0];

  std::complex<T> coeffs[2] = {
    (eigens[1] * yo - dyo) / denom,
    (-eigens[0] * yo + dyo) / denom
  };

  // Stub for determining particular solution
  // TODO: Find particular solution
  (void) f;
  endomorphism<T> particular = zeroFunction<T>;

  return endomorphism<T>(
    /**
     * @param t Temporal input
     * @return Position of mass in system at t
     */
    [=] (T t) {
      /*
       * Since the characteristic equation has real coefficients, both
       * eigenvalues will be complex conjugates of each other. Since
       * conjugation is preserved under exponentiation, we are adding
       * two conjugates and thus cancelling the imaginary part. This
       * means that we only have to take the real part of the solution
       * (assuming that the particular solution has real range).
       */
      return (
        coeffs[0] * std::exp(eigens[0] * t) +
        coeffs[1] * std::exp(eigens[1] * t) +
        particular(t)
      ).real();
    }
  );
}
{% endhighlight %}
