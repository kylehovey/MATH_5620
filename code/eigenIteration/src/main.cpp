#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  auto A = Mtx::hilbert(5);

  const auto [ bigEigen, x ] = A.largestEigenpair();
  const auto [ smolEigen, y ] = A.smallestEigenpair();
  const auto condition = A.conditionNumber();

  std::cout << "A\n";
  std::cout << A << std::endl;
  std::cout << "x vector\n";
  std::cout << x << std::endl;
  std::cout << "y vector\n";
  std::cout << y << std::endl;
  std::cout << "Condition Number of A:" << std::endl;
  std::cout << condition << std::endl;

  std::cout << "Largest Eigenvalue\n";
  std::cout << bigEigen << std::endl;
  std::cout << std::endl;
  std::cout << "A * x\n";
  std::cout << A * x << std::endl;
  std::cout << "Largest Eigenvalue * x\n";
  std::cout << bigEigen * x << std::endl;
  std::cout << "Smallest Eigenvalue\n";
  std::cout << smolEigen << std::endl;
  std::cout << std::endl;
  std::cout << "A * y\n";
  std::cout << A * y << std::endl;
  std::cout << "Smallest Eigenvalue * y\n";
  std::cout << smolEigen * y << std::endl;

  return EXIT_SUCCESS;
}
