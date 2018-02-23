#include "../../matrix/src/matrix/matrix.h"
#include <functional>
#include <vector>
#include <iostream>
#include <limits>

using Mtx = Matrix::Matrix<double>;

template <typename T>
using dubs = std::tuple<std::tuple<T, T>, std::tuple<T, T>>;

template <typename T>
using planeToScalar = std::function<T(T, T)>;

template <typename T>
using stencilGen = std::function<
  std::vector<
    std::tuple<
      T,
      std::tuple<T, T>
    >
  >
>;

/**
 * Solve ∇²u(x, y) = f(x, y)
 * @param size Mesh size
 * @param domain Lower-left and upper-right coordinates of domain
 * @param stencil Function that takes tuple of evaluation point coordinates
 *  and returns a vector of tuples for the resulting stencil that also
 *  includes multipliers. Example:
 *    { (-4, (2, 2)), (1, (2, 3)), (1, (2, 1)) }
 * @param dirichlet Function u(x, y) along boundary
 * @return Matrix of values that resemble a solution over u's domain
 */
template <typename T>
Matrix::Matrix<T> solveLaplace(
    const uint& size,
    const dubs<T>& domain,
    const stencilGen<T>& stencil,
    const planeToScalar<T>& dirichlet
) {
  (void) size;
  (void) domain;
  (void) stencil;
  (void) dirichlet;
}

int main() {
  // Define domain
  const dubs<double> domain = { { 0, 0 }, { 1, 1 } };

  // G(x, y) along the boundaries
  const planeToScalar<double> dOmega =
    [](const double& a, const double& b) -> double {
      (void) b;

      return a == 1.0 ? 5.0 : 0;
    };

  return EXIT_SUCCESS;
}
