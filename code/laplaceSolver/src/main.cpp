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
using coord = std::tuple<T, T>;

template <typename T>
using stencil = std::vector<std::tuple<T, coord<T>>>;

template <typename T>
using stencilGen = std::function<stencil<T>(
    const std::tuple<T, T>&,
    const T& h
)>;

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
    const stencilGen<T>& makeStencil,
    const planeToScalar<T>& dirichlet
) {
  (void) size;
  (void) domain;
  (void) makeStencil;
  (void) dirichlet;

  return Matrix::Matrix<T>(size, size);
}

int main() {
  // Define domain
  const dubs<double> domain = { { 0, 0 }, { 1, 1 } };

  // G(x, y) along the boundaries
  const planeToScalar<double> boundary =
    [](const double& a, const double& b) -> double {
      (void) b;

      return a == 1.0 ? 5.0 : 0;
    };

  // Define size of mesh
  const uint size = 5;

  // Stencil generation
  const stencilGen<double> fivePoint = [](
      const coord<double>& center,
      const double& h
  ) {
    const auto [ x, y ] = center;

    return stencil<double>({
      { -4,{ x + 0, y + 0 } },
      { 1, { x + h, y + 0 } },
      { 1, { x - h, y + 0 } },
      { 1, { x + 0, y + h } },
      { 1, { x + 0, y - h } }
    });
  };

  const stencilGen<double> ninePoint = [](
      const coord<double>& center,
      const double& h
  ) {
    const auto [ x, y ] = center;

    return stencil<double>({
      { -8, { x - 0, y + 0 } },
      { 1,  { x - 0, y + h } },
      { 1,  { x - 0, y - h } },
      { 1,  { x + h, y + 0 } },
      { 1,  { x + h, y + h } },
      { 1,  { x + h, y - h } },
      { 1,  { x - h, y + 0 } },
      { 1,  { x - h, y + h } },
      { 1,  { x - h, y - h } }
    });
  };

  const auto soln = solveLaplace<double>(size, domain, fivePoint, boundary);

  std::cout << soln << std::endl;

  return EXIT_SUCCESS;
}
