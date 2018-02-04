#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx A({
      { 0, -2, 3, 4 },
      { 2, 4, 1, 9 },
      { -3, 2, 0, 2 },
      { 3, 2, 1, 2 }
  });

  const auto [ P, L, U ] = A.LUFactorize();

  const auto res = P * L * U;

  std::cout << A << std::endl;
  std::cout << P << std::endl;
  std::cout << L << std::endl;
  std::cout << U << std::endl;
  std::cout << res << std::endl;

  return EXIT_SUCCESS;
}
