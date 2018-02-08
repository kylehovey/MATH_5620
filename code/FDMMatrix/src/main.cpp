#include "../../matrix/src/matrix/matrix.h"
#include <iostream>
#include <vector>
#include <cmath>

using Mtx = Matrix::Matrix<double>;

/**
 * Solve u'' = f given boundary conditions
 * @param a Left boundary
 * @param b Right boundary
 * @param ua u(a)
 * @param ub u(b)
 * @param f Driving function
 * @param n Size of mesh
 * @return Column vector of solution
 */
template <typename T>
Matrix::Matrix<T> solveElliptic(
    const T& a,
    const T& b,
    const T& ua,
    const T& ub,
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
  const auto D = Matrix::Matrix<T>::genFDMatrix(n, 2);

  // Solve for solution vector
  const auto soln = Matrix::Matrix<T>::solve(
      D,
      F - init,
      Matrix::Solve::Thompson
   );

  return soln;
}

int main() {
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

  // Generate exact solution vector
  std::function<double(double)> exact = [](double x) {
    return 2.5 * x + 2.5 - std::sin(M_PI * x) / (M_PI * M_PI);
  };

  const double h = (b - a) / (double) meshSize;
  Mtx uExact(meshSize, 1, [&](const uint& i, const uint& j) {
      (void) j;
      return exact(a + i * h);
  });

  const auto soln = solveElliptic<double>(a, b, ua, ub, f, meshSize);

  std::cout << "Solved solution\n";
  std::cout << soln << std::endl;

  std::cout << "Exact solution\n";
  std::cout << uExact << std::endl;

  std::cout << "Error vector\n";
  const auto E = soln - uExact;
  std::cout << E << std::endl;

  std::cout << "1-norm of error vector\n";
  std::cout << Matrix::Matrix<double>::vNorm(E, 1) << std::endl;

  std::cout << "2-norm of error vector\n";
  std::cout << Matrix::Matrix<double>::vNorm(E, 2) << std::endl;

  std::cout << "infinity-norm of error vector\n";
  std::cout << Matrix::Matrix<double>::vNorm(E, std::numeric_limits<uint>::max()) << std::endl;

  return EXIT_SUCCESS;
}
