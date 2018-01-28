#include "matrix/matrix.h"
#include <iostream>
#include <vector>

int main() {
  Matrix::Matrix<double> A({
      { 1, 0, 0, 0, 0 },
      { 0, 1, 0, 0, 0 },
      { 0, 0, 1, 0, 0 },
      { 0, 0, 0, 1, 0 },
      { 0, 0, 0, 0, 1 }
  });

  Matrix::Matrix<double> B({
      { 1, 2, 3, 4, 5 },
      { 6, 7, 8, 9, 0 },
      { 1, 2, 3, 4, 5 },
      { 6, 7, 8, 9, 0 },
      { 6, 7, 8, 9, 0 }
  });

  std::cout << A.isDiagonal() << std::endl;

  auto C = 5 * (A + B) * A;

  if (C != A) {
    std::cout << "Not equal!" << std::endl;
  }

  std::cout << A.trace() << std::endl;

  auto D = C;

  if (C == D) {
    std::cout << "Equal!" << std::endl;
  }

  std::cout << A << std::endl;
  std::cout << B << std::endl;
  std::cout << C << std::endl;

  return EXIT_SUCCESS;
}
