#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  for (uint i = 0; i < 10; ++i) {
    const auto size = (i + 1) * 10;
    auto A = Mtx::genFDMatrix(size, 2);
    const auto smallEigen = std::get<0>(A.smallestEigenpair());
    const auto condition = A.conditionNumber();

    std::cout << "2nd order FDM matrix smallest eigenvalue for size ";
    std::cout << size << ":" << std::endl;
    std::cout << smallEigen << std::endl;

    std::cout << "2nd order FDM matrix condition number for size ";
    std::cout << size << ":" << std::endl;
    std::cout << condition << std::endl;
  }

  return EXIT_SUCCESS;
}
