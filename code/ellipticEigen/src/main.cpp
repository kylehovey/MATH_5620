#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  for (uint i = 0; i < 15; ++i) {
    const auto size = (i + 1) * 5;
    auto A = Mtx::genFDMatrix(size, 2);
    const auto bigEiegen = std::get<0>(A.largestEigenpair());
    std::cout << "2nd Order FDM Matrix Spectral Radius for size ";
    std::cout << size << ":" << std::endl;
    std::cout << bigEiegen << std::endl;
  }

  return EXIT_SUCCESS;
}
