#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <limits>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx v({{1},{2},{3},{4}});

  auto oneNorm = Mtx::vNorm(v, 1);
  auto twoNorm = Mtx::vNorm(v, 2);
  auto infNorm = Mtx::vNorm(v, std::numeric_limits<uint>::max());

  std::cout << "1-norm of vector:" << std::endl;
  std::cout << oneNorm << std::endl;
  std::cout << "2-norm of vector:" << std::endl;
  std::cout << twoNorm << std::endl;
  std::cout << "Infinity-norm of vector:" << std::endl;
  std::cout << infNorm << std::endl;

  return EXIT_SUCCESS;
}
