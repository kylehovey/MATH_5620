#include "matrix/matrix.h"
#include <iostream>
#include <vector>

int main() {
  Matrix::Matrix<int> A(5, 5, (Matrix::binaryDual<int>) [](const uint& a, const uint& b) {
    return a % (b + 1);
  });

  Matrix::Matrix<int> B({
      { 1, 2, 3, 4, 5 },
      { 6, 7, 8, 9, 0 },
      { 1, 2, 3, 4, 5 },
      { 6, 7, 8, 9, 0 },
      { 6, 7, 8, 9, 0 }
  });

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << A.getVal(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << B.getVal(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << A.getVal(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << A.trace() << std::endl;


  return EXIT_SUCCESS;
}
