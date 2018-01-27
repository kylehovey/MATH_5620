#include "matrix/matrix.h"
#include <iostream>

int main() {
  Matrix::Matrix<double> A(5, 5, (Matrix::binaryDual<int>) [](const uint& a, const uint& b) {
      return a - b;
  });

  Matrix::Matrix<int> B(5, 5);

  //for (int i = 0; i < 5; i++) {
  //  for (int j = 0; j < 5; j++) {
  //    std::cout << A.getVal(i, j) << " ";
  //  }
  //  std::cout << std::endl;
  //}

  A.transpose();

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << B.getVal(i, j) << " ";
    }
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
