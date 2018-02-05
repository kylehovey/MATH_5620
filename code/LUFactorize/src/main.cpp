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

  Mtx _x({{1},{2},{3},{4}});

  auto b = A * _x;

  auto x = Mtx::solve(A, b, Matrix::Solve::LU);

  std::cout << x << std::endl;

  return EXIT_SUCCESS;
}
