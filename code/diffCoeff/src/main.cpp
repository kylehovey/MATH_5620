#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  std::cout << "1st derivative (m = 1) with order 4 accuracy (n = 4)\n";
  const auto coeffs = Mtx::genFDCoeff(2, 2);

  for (auto i = 0u; i < coeffs.size(); ++i) {
    std::cout << coeffs[i] << std::endl;
  }

  return EXIT_SUCCESS;
}
