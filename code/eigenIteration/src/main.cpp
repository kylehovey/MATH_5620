#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  auto A = Mtx::hilbert(5);

  const auto [ bigEigen, x ] = A.largestEigenpair();

  std::cout << "A\n";
  std::cout << A << std::endl;

  std::cout << "Largest Eigenvalue\n";
  std::cout << bigEigen << std::endl;
  std::cout << std::endl;
  std::cout << "x vector\n";
  std::cout << x << std::endl;
  std::cout << "A * x\n";
  std::cout << A * x << std::endl;
  std::cout << "Largest Eigenvalue * x\n";
  std::cout << bigEigen * x << std::endl;

  return EXIT_SUCCESS;
}
