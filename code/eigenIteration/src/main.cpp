#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx A({
      { -2, 1, 0, 0 },
      { 1, -2, 1, 0 },
      { 0, 1, -2, 1 },
      { 0, 0, 1, -2 }
  });

  Mtx x({{1},{2},{3},{4}});

  for (auto i = 0u; i < 100; ++i) {
    auto b = A * x;
    double mult = 1.0 / Mtx::vNorm(b, 2);
    x = mult * b;
  }

  std::cout << "A\n";
  std::cout << A << std::endl;
  std::cout << "x vector\n";
  std::cout << x << std::endl;

  auto axe = A * x;
  auto eigenVal = axe.getVal(0, 0) / x.getVal(0, 0);

  std::cout << "Eigenvalue\n";
  std::cout << eigenVal << std::endl;
  std::cout << std::endl;
  std::cout << "A * x\n";
  std::cout << A * x << std::endl;
  std::cout << "Eigenvalue * x\n";
  std::cout << eigenVal * x << std::endl;

  return EXIT_SUCCESS;
}
