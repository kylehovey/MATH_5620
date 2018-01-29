#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>

using Mtx = Matrix::Matrix<double>;

int main() {
  Mtx A({ { 8, -2, 3 }, { 2, 5, 1 }, { -3, 2, 6 } });

  Mtx invD(3, 3, [&](const uint& a, const uint& b) {
    return a == b ? 1 / A.getVal(a, b) : 0;
  });

  auto R = A.lTriangular() + A.uTriangular();

  Mtx b({ { 7 }, { 4 }, { -10 } });

  auto x = b;

  for (unsigned int i = 0; i < 500; ++i) {
    x = invD * (b - R * x);
  }

  std::cout << x << std::endl;

  return EXIT_SUCCESS;
}
