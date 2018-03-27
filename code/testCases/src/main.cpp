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
