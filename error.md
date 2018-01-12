---
math: true
permalink: /error
layout: default
---

# Computing Hardware Epsilon

**Routine Name**: Relative and Absolute Error

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This templated function is designed to calculate absolute and relative error when the exact form of a number is known. If we define \\(v\\) to be the exact amount and \\(v_{\text{approx}}\\) to be an approximation, then absolute error (\\( \epsilon \\)) is defined to be:

\\[ \epsilon := \Big \lvert v - v_{\text{approx}} \Big \rvert \\]

and relative error is defined to be:

\\[ \eta := \frac{\epsilon}{\lvert v \rvert} \\]

**Input**:

The first parameter will be the approximation, and the second parameter will be the exact value. A type is also specified in the template function.

**Output**:

Either the absolute or relative error, depending on the function call.

**Usage/Example**:

To use this function, you have to include it in the file you need to call it from.

{% highlight C++ %}
#include "absolute_error/absolute_error.h"
#include "relative_error/relative_error.h"
#include <limits>
#include <iostream>

int main() {
  const double approx = 3.2;
  const double exact = M_PI;

  std::cout << "PI: " << exact << '\n';
  std::cout << "Approximation: " << approx << '\n';
  std::cout << "- - - - - - - - - -" << '\n';
  std::cout << "Absolute: " << absoluteError<double>(approx, exact) << '\n';
  std::cout << "Relative: " << relativeError<double>(approx, exact) << '\n';

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

{% highlight C++ %}
PI: 3.14159
Approximation: 3.2
- - - - - - - - - -
Absolute: 0.0584073
Relative: 0.0185916
{% endhighlight %}

**Implementation/Code:**

Computing absolute error is straightforward:

{% highlight c++ %}
template <typename T>
/**
 * @param approximate The value that is approximated
 * @param exact The exact value to be compared against
 * @return The relative error
 */
T relativeError(const T approximate, const T exact) {
  return absoluteError<T>(approximate, exact) / exact;
}
{% endhighlight %}

Computing relative error utilizes the computation for absolute error:

{% highlight c++ %}
template <typename T>
/**
 * @param approximate The value that is approximated
 * @param exact The exact value to be compared against
 * @return The relative error
 */
T relativeError(const T approximate, const T exact) {
  return absoluteError<T>(approximate, exact) / exact;
}
{% endhighlight %}
