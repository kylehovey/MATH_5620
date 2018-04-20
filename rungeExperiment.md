---
math: true
permalink: /rungeExperiment
title: Runge Experiment
layout: page
---

**Routine Name**: Runge Experiment

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This code tests the accuracy of a given method against some examples found in the course textbook. In \\(7.1\\), we analyze \\(u' = -\sin(x); u_0 = 1; dt = 10^{-3}\\). In \\(7.2\\), we change our driving function to \\(u' = \lambda (u - \cos(x)) - \sin(x)\\) and require \\(\lambda = -10\\). And finally, for \\(7.3\\), we require \\(\lambda = -2100\\).

**Input**: A method for solving the examples

**Output**: The solutions given by the method

**Usage/Example**:

{% highlight C++ %}
#include <iostream>
#include <cmath>
#include "../../rungeKutta/src/rungeKutta/rungeKutta.h"

int main() {
  /*
   * 7.1
   */
  const auto dt = 1E-3;
  const auto uo = 1;
  const auto simple = [](const double& t, const double& u) -> double {
    (void) u;
    return -std::sin(t);
  };
  const auto T = 2;

  auto soln = RungeKutta::genOrderFourSolution<double>(simple, dt, uo);

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

  soln = RungeKutta::genOrderFourSolution<double>(modified, dt, uo);

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

  soln = RungeKutta::genOrderFourSolution<double>(largeLambda, dt, uo);

  std::cout << "7.3:" << std::endl;
  std::cout << "U_2000 = approx(2): " << soln(T) << std::endl;
  std::cout << "Exact(2): " << std::cos(T) << std::endl;
  std::cout << "Error: " << std::abs(soln(T) - std::cos(T)) << std::endl;

  return EXIT_SUCCESS;
}
{% endhighlight %}

Output:

The output for the first two examples *still* has issues. But it is still easy to see that the error blows up for \\(7.3\\) where \\(\lambda = -21000\\). In fact, it blows up so much that it overflows and results in a `nan` result.

{% highlight C++ %}
7.1:
U_2000 = approx(2): -0.818512
Exact(2): -0.416147
Error: 0.402365

7.2:
U_2000 = approx(2): -0.507163
Exact(2): -0.416147
Error: 0.0910165

7.3:
U_2000 = approx(2): -0.416682
Exact(2): -0.416147
Error: 0.000534717
{% endhighlight %}

Decreasing the step size to \\(10^{-4}\\) yields accurate results for the example in \\(7.3\\). This shows that the stability region for a fourth-order Runge-Kutta method also has a tight tolerance around \\(10^{-3}\\) as a time step.

{% highlight C++ %}
7.1:
U_2000 = approx(2): -0.818512
Exact(2): -0.416147
Error: 0.402365

7.2:
U_2000 = approx(2): -0.507163
Exact(2): -0.416147
Error: 0.0910165

7.3:
U_2000 = approx(2): -0.416682
Exact(2): -0.416147
Error: 0.000534717
{% endhighlight %}

**Implementation/Code:**

All solver implementation for this assignment was covered in the last assignment.
