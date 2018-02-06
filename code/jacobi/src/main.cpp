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

  Mtx _x({{1},{2},{3},{4}});

  auto b = A * _x;

  auto x = Mtx::solve(A, b, Matrix::Solve::Jacobi);

  std::cout << x << std::endl;

  return EXIT_SUCCESS;
}
