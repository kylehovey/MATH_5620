---
math: true
permalink: /testCases
title: Some Test Cases
layout: page
---

**Routine Name**: Exponential Lambda Solution

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This method generates the analytic solution to the basic differential equation:

# \\[ u' = \lambda u; u(0) \equiv \alpha \\]

which has the solution

# \\[ u(t) = \alpha e^{\lambda t} \\]

**Input**:

* \\( \alpha \\)
* \\( \lambda \\)

**Output**:

A lambda function that can be evaluated for an exact solution to the analytic equation (up to roundoff error).

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include "testCases/testCases.h"

int main() {
  const auto alpha = 10.0;
  const auto lambda = -1.0;

  const auto soln = TestCases::genLambdaSolution(lambda, alpha);

  for (auto i = 0u; i < 10; ++i) {
    std::cout << "u(" << i << ") = " << soln(i) << std::endl;
  }

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
u(0) = 10
u(1) = 3.67879
u(2) = 1.35335
u(3) = 0.497871
u(4) = 0.183156
u(5) = 0.0673795
u(6) = 0.0247875
u(7) = 0.00911882
u(8) = 0.00335463
u(9) = 0.0012341
{% endhighlight %}

**Implementation/Code:**

{% highlight c++ %}
#include <functional>
#include <cmath>

template <typename T>
using endomorphism = std::function<T(const T&)>;

namespace TestCases {
  /**
   * Generates the analytic solution to the differential equation
   *  u' = λu; u(0) := α
   * @param lambda Lambda in the equation above
   * @param alpha Initial condition at u(0)
   * @return Analytic solution
   */
  template <typename T>
  endomorphism<T> genLambdaSolution(const T& lambda, const T& alpha) {
    return [=](const T& t) -> T {
      return alpha * std::exp(lambda * t);
    };
  }
};
{% endhighlight %}

<hr>

**Routine Name**: Logistic Equation (From Assignment 1)

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

A logistic equation is a solution \\(P\\) to the differential equation:

## \\[ \frac{dP}{dt} = \gamma P - \beta P^2 \\]

The solution of which can be massaged to the form:

## \\[ \frac{\gamma}{(\frac{\gamma}{P_o} - \beta)e^{-\gamma t} + \beta} \ni P_o := P(0)\\]

This code takes the free parameters for the logistic equation and returns a function that takes a templated time variable and returns a value of the same type at time `t`.

**Input**:

The free parameters that define the logistic equation: \\( \\{ \gamma, \beta, P_o \\} \\), as well as a type `T` that the function should be defined over.

**Output**:

A function that takes a parameter `t` of type `T` and returns the value of the specified logisic equation evaluated at `t`.

**Usage/Example**:

To use this function, you have to include it in the file you need to call it from.

{% highlight C++ %}
#include "logistic/logistic.h"
#include <limits>
#include <iostream>

int main() {
  // Parameters for logistic equation
  const double
    gamma = 2,
    beta = 1,
    Po = 1;

  // Generate the required function
  auto logistic = genLogistic<double>(gamma, beta, Po);

  // Call it for some basic values
  for (int i = -10; i <= 10; i++) {
    std::cout << i << " -> " << logistic(i) << std::endl;
  }

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
-10 -> 4.12231e-09
-9 -> 3.046e-08
-8 -> 2.2507e-07
-7 -> 1.66306e-06
-6 -> 1.22883e-05
-5 -> 9.07957e-05
-4 -> 0.0006707
-3 -> 0.00494525
-2 -> 0.0359724
-1 -> 0.238406
0 -> 1
1 -> 1.76159
2 -> 1.96403
3 -> 1.99505
4 -> 1.99933
5 -> 1.99991
6 -> 1.99999
7 -> 2
8 -> 2
9 -> 2
10 -> 2
{% endhighlight %}

**Implementation/Code:**

I made this function as abstract as possible. I probably could have gone with partial function application to bind \\( \gamma \\), \\( \beta \\), and \\( P_o \\), but C++ does not have a very clean way of doing this yet. Instead, I created a small factory for the logistic equation that took the free parameters and returns a function that is a logistic map.

{% highlight c++ %}
template <typename T>
using endomorphism = std::function<T(T)>;

template <typename T>
endomorphism<T> genLogistic(T gamma, T beta, T Po) {
  auto amp = ( gamma / Po - beta);

  return endomorphism<T>(
    [=] (T t) {
      return gamma / (amp * exp(-gamma * t) + beta);
    }
  );
}
{% endhighlight %}
