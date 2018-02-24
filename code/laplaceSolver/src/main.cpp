#include "../../matrix/src/matrix/matrix.h"
#include <functional>
#include <vector>
#include <iostream>
#include <limits>

using Mtx = Matrix::Matrix<double>;

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
    const coord<T>& domain,
    const stencilGen<T>& makeStencil,
    const planeToScalar<T>& dirichlet
) {
  (void) domain;
  (void) makeStencil;
  (void) dirichlet;

  // Find limits of domain
  const auto a = std::get<0>(domain);
  const auto b = std::get<1>(domain);
  const auto h = (b - a) / (size - 1);

  // Internal limits of sampling
  const auto _a = a + h;
  const auto _b = a + b - h;

  // Instantiate output grid (interior of mesh)
  Matrix::Matrix<T> rhs(
      size - 2,
      size - 2,
      [&](const uint& row, const uint& col) {
        // Determine location in domain
        const auto x = a + h * (col + 1);
        const auto y = a + h * (row + 1);

        // Initialize an accumulator
        T acc = 0;

        // Determine stencil
        const auto samples = makeStencil({ x, y }, h);

        // For each point in stencil, subtract off boundary
        for (const auto& sample : samples) {
          const auto [ mult, coords ] = sample;
          const auto [ _x, _y ] = coords;

          // If outside of interior
          if (_x < _a || x > _b || _y < _a || _y > _b) {
            acc -= dirichlet(_x, _y);
          }
        }

        return acc;
      }
  );

  return rhs;
}

int main() {
  // Define domain
  const coord<double> domain = { 0, 1 };

  // G(x, y) along the boundaries
  const planeToScalar<double> boundary =
    [](const double& x, const double& y) -> double {
      (void) x;

      return y == 1 ? 5.0 : 0;
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
