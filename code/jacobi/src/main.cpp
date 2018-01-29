#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

int main() {
  // Input matrix
  Matrix::Matrix<double> A({
      { 8, -2, 3 },
      { 2, 5, 1 },
      { -3, 2, 6 }
  });

  // D and R matrices
  Matrix::Matrix<double> invD(3, 3, [&](const uint& a, const uint& b) {
    return a == b ? 1.0 / A.getVal(a, b) : 0;
  });

  Matrix::Matrix<double> R(3, 3, [&](const uint& a, const uint& b) {
      return a != b ? A.getVal(a, b) : 0;
  });

  // b vector
  Matrix::Matrix<double> b({
      { 7 },
      { 4 },
      { -10 }
  });

  // x vector
  Matrix::Matrix<double> x = b;

  for (unsigned int i = 0; i < 500; ++i) {
    x = invD * (b - (R * x));
  }

  std::cout << x << std::endl;

  return EXIT_SUCCESS;
}
