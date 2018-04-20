---
math: true
permalink: /implicitExperiment
title: Implicit Experiment
layout: page
---

**Routine Name**: Implicit Experiment

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This code goes over the exact same model problems discussed in problem one on this assignment, but now it uses an implicit setup to evaluate the examples.

**Input**: A method for solving the examples

**Output**: The solutions given by the method

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include <cmath>
#include "../../euler/src/euler/euler.h"

int main() {
  /*
   * 7.1
   */
  const auto dt = 1E-6;
  const auto uo = 1;
  const auto simple = [](const double& t, const double& u) -> double {
    (void) u;
    return -std::sin(t);
  };
  const auto T = 2;

  auto soln = Euler::genExplicitEulerSolution<double>(simple, dt, uo);

  std::cout << "7.1:" << std::endl;
  std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
  std::cout << "Exact(2): " << std::cos(T) << std::endl;
  std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;
  std::cout << std::endl;

  /*
   * 7.2
   */
  const auto modified = [](const double& t, const double& u) -> double {
    return -10 * (u - std::cos(t)) - std::sin(t);
  };

  soln = Euler::genExplicitEulerSolution<double>(modified, dt, uo);

  std::cout << "7.2:" << std::endl;
  std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
  std::cout << "Exact(2): " << std::cos(T) << std::endl;
  std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;
  std::cout << std::endl;

  /*
   * 7.3
   */
  const auto largeLambda = [](const double& t, const double& u) -> double {
    return -2100 * (u - std::cos(t)) - std::sin(t);
  };

  soln = Euler::genExplicitEulerSolution<double>(largeLambda, dt, uo);

  std::cout << "7.3:" << std::endl;
  std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
  std::cout << "Exact(2): " << std::cos(T) << std::endl;
  std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

Notice how the output for the first two examples still has very high error. I'm still not completely sure why this happened (as there was very little error when testing the exact same code for the last assignment), but it is easy to see that the error is very small for \\(7.3\\), which is very different from the case with Explicit Euler. This is because Implicit Euler is absolutely stable for a much more comprehensive region than Explicit Euler.

{% highlight C++ %}
7.1:
U_2000 = approx(2): -0.818595
Exact(2): -0.416147
Error: 0.402448

7.2:
U_2000 = approx(2): -0.506977
Exact(2): -0.416147
Error: 0.0908299

7.3:
U_2000 = approx(2): -0.41658
Exact(2): -0.416147
Error: 0.000432812
{% endhighlight %}

We can also increase the step size and see that Implicit Euler still remains convergent as it is much more stable than Explicit Euler.

{% highlight C++ %}
USING STEP SIZE: 0.01
==================

7.1:
U_2000 = approx(2): -0.818595
Exact(2): -0.416147
Error: 0.402448

7.2:
U_2000 = approx(2): -0.506977
Exact(2): -0.416147
Error: 0.0908303

7.3:
U_2000 = approx(2): -0.41658
Exact(2): -0.416147
Error: 0.000432724

USING STEP SIZE: 0.1
==================

7.1:
U_2000 = approx(2): -0.818595
Exact(2): -0.416147
Error: 0.402448

7.2:
U_2000 = approx(2): -0.506985
Exact(2): -0.416147
Error: 0.0908378

7.3:
U_2000 = approx(2): -0.41658
Exact(2): -0.416147
Error: 0.000432848
{% endhighlight %}

**Implementation/Code:**

All solver implementation for this assignment was covered in the last assignment.
