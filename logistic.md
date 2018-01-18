---
math: true
permalink: /logistic
title: Logistic Equation Analytic Solution
layout: page
---

**Routine Name**: Logistic Equation

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

A logistic equation is a solution \\(P\\) to the differential equation:

## \\[ \frac{dP}{dt} = \alpha P - \beta P^2 \\]

The solution of which can be massaged to the form:

## \\[ \frac{\alpha}{(\frac{\alpha}{P_o} - \beta)e^{-\alpha t} + \beta} \ni P_o := P(0)\\]

This code takes the free parameters for the logistic equation and returns a function that takes a templated time variable and returns a value of the same type at time `t`.

**Input**:

The free parameters that define the logistic equation: \\( \\{ \alpha, \beta, P_o \\} \\), as well as a type `T` that the function should be defined over.

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
    alpha = 2,
    beta = 1,
    Po = 1;

  // Generate the required function
  auto logistic = genLogistic<double>(alpha, beta, Po);

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

I made this function as abstract as possible. I probably could have gone with partial function application to bind \\( \alpha \\), \\( \beta \\), and \\( P_o \\), but C++ does not have a very clean way of doing this yet. Instead, I created a small factory for the logistic equation that took the free parameters and returns a function that is a logistic map.

{% highlight c++ %}
template <typename T>
using homomorphism = std::function<T(T)>;

template <typename T>
homomorphism<T> genLogistic(T alpha, T beta, T Po) {
  auto amp = ( alpha / Po - beta);

  return homomorphism<T>(
    [=] (T t) {
      return alpha / (amp * exp(-alpha * t) + beta);
    }
  );
}
{% endhighlight %}
