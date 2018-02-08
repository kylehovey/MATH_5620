#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  std::cout << "2nd derivative (m = 2) with order 2 accuracy (n = 4)\n";

  const auto coeffs = Mtx::genFDCoeff(2, 4);

  for (auto i = 0u; i < coeffs.size(); ++i) {
    std::cout << coeffs[i] << std::endl;
  }

  return EXIT_SUCCESS;
}
