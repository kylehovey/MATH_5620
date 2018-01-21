---
math: true
permalink: /epsilon
title: Hardware Epsilon
layout: page
---

**Routine Name**: Hardware Epsilon

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This templated function is designed to compute the hardware epsilon on any machine using an abstract type defined per-use-case.

From [Wikipedia](https://en.wikipedia.org/wiki/Machine_epsilon):

>_Machine epsilon gives an upper bound on the relative error due to rounding in floating point arithmetic. This value characterizes computer arithmetic in the field of numerical analysis, and by extension in the subject of computational science. The quantity is also called macheps or unit roundoff, and it has the symbols Greek epsilon \\( \epsilon \\) or bold Roman **\\( u \\)**, respectively._

**Input**:

Either nothing, or a starting value to compute epsilon from.

**Output**:

The machine epsilon for the computer that the code runs on.

**Usage/Example**:

To use this function, you have to include it in the file you need to call it from.

{% highlight C++ %}
#include "epsilon/epsilon.h"
#include <limits>
#include <iostream>

int main() {
  std::cout << computeEpsilon<double>() << std::endl;
  std::cout << std::numeric_limits<double>::epsilon() << std::endl;

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
2.22045e-16
2.22045e-16
{% endhighlight %}

**Implementation/Code:**

To find the smallest quanta represented in a floating type `T`, I wrote a templated function that will reduce the exponent (by dividing by two) of an `epsilon` value of type `T` until it becomes an identity element under floating point addition. To ensure that it is the identity, I will add it to an integer and ensure that it is unchanged after the operation:

{% highlight c++ %}
template <typename T>
/**
 * @param eps Initial value for epsilon
 * @return Machine epsilon fir given type
 */
T computeEpsilon(T eps = 1) {
  // While epsilon is still visible to addition
  while (1 + eps != 1) {
    // Reduce exponent by one
    eps /= 2;
  }

  // Eps is no longer visible, show it's value before it was halved
  return eps * 2;
}
{% endhighlight %}
