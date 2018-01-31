#include "matrix/matrix.h"
#include <iostream>
#include <vector>

int main() {
  auto A = Matrix::Matrix<int>(5, 5, [](const uint& a, const uint& b) {
    return a - b;
  });

  std::cout << A << std::endl;

  return EXIT_SUCCESS;
}
