#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using Mtx = Matrix::Matrix<double>;

/**
 * Solve u'' = f given boundary conditions
 * @param a Left boundary
 * @param b Right boundary
 * @param ua u(a)
 * @param ub u(b)
 * @param k Vector representing k function
 * @param f Driving function
 * @param n Size of mesh
 * @return Column vector of solution
 */
template <typename T>
Matrix::Matrix<T> solveEllipticWithK(
    const T& a,
    const T& b,
    const T& ua,
    const T& ub,
    const Matrix::Matrix<T>& k,
    const std::function<T(T)>& f,
    const uint& n
) {
  (void) ua;
  (void) ub;

  // Generate F vector
  const T h = (b - a) / (T) n;
  const Matrix::Matrix<T> F(n, 1, [&](const uint& i, const uint& j) -> T {
      (void) j;
      return std::pow(h, 2) * f(a + i * h);
  });

  // Use initial conditions
  Matrix::Matrix<T> init(n, 1);
  init.setVal(0, 0, ua);
  init.setVal(n - 1, 0, ub);

  // Generate differential operator for second derivative
  auto D = Matrix::Matrix<T>::genFDMatrix(n, 2);

  // Modify it so that it includes the derivative for k(x)
  D.fillWith([&](const uint& a, const uint& b) -> T {
      const auto val = D.getVal(a, b);
      return val * k.getVal(b, 0);
  });

  // Solve for solution vector
  const auto soln = Matrix::Matrix<T>::solve(
      D,
      F - init,
      Matrix::Solve::Thompson
   );

  return soln;
}

int main() {
  // Perfectly secure seed
  srand(time(NULL));

  // Mesh size
  const uint meshSize = 10;

  // Limits
  const double a = 0;
  const double b = 1;

  // Function u at limits
  const double ua = 2.5;
  const double ub = 5.0;

  // Driving function
  std::function<double(double)> f = [](double x) {
    return std::sin(M_PI * x);
  };

  // K matrix
  const Mtx k(meshSize, 1, [](const uint& a, const uint& b) {
      (void) a;
      (void) b;
      return rand() % 40 + 10;
  });

  const auto soln = solveEllipticWithK<double>(a, b, ua, ub, k, f, meshSize);

  std::cout << "Solved solution\n";
  std::cout << soln << std::endl;

  return EXIT_SUCCESS;
}
