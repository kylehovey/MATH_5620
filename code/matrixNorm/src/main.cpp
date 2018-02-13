#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <limits>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx A({
      {1, 5, 1, 0},
      {2, 0, 0, 1},
      {3, 8, 5, 2},
      {4, 9, 2, 1}
  });

  auto oneNorm = Mtx::mNorm(A, 1);
  auto infNorm = Mtx::mNorm(A, std::numeric_limits<uint>::max());

  std::cout << "1-norm of matrix:" << std::endl;
  std::cout << oneNorm << std::endl;
  std::cout << "Infinity-norm of matrix:" << std::endl;
  std::cout << infNorm << std::endl;

  return EXIT_SUCCESS;
}
