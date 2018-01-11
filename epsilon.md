---
math: true
permalink: /epsilon
layout: default
---

# Computing Hardware Epsilon

From Wikipedia:

>_Machine epsilon gives an upper bound on the relative error due to rounding in floating point arithmetic. This value characterizes computer arithmetic in the field of numerical analysis, and by extension in the subject of computational science. The quantity is also called macheps or unit roundoff, and it has the symbols Greek epsilon \\( \epsilon \\) or bold Roman **\\( u \\)**, respectively._

To find the smallest quanta represented in a floating type `T`, I wrote a templated function that will reduce the exponent (by dividing by two) of an `epsilon` value of type `T` until it becomes an identity element under floating point addition. To ensure that it is the identity, I will add it to an integer and ensure that it is unchanged after the operation:

{% highlight c++ %}
template <typename T>
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

Using this code, we can compare it against the built-in numeric limits utility function for determining hardware epsilon in the C++ standard library.

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

>2.22045e-16 <br />
>2.22045e-16

This verifies the machine epsilon computation. \\( \blacksquare \\)
