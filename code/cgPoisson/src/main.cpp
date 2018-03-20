#include "../../matrix/src/matrix/matrix.h"
#include "../../image/src/image/image.h"
#include <functional>
#include <vector>
#include <iostream>
#include <limits>
#include <cmath>

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
 * @param driver Driving function (f(x, y) on the r.h.s.)
 * @param stencil Function that takes tuple of evaluation point coordinates
 *  and returns a vector of tuples for the resulting stencil that also
 *  includes multipliers. Offsets should be integer multiples of h. Example:
 *    { (-4, (2, 2)), (1, (2, 3)), (1, (2, 1)) }
 * @param dirichlet Function u(x, y) along boundary
 * @return Matrix of values that resemble a solution over u's domain
 */
template <typename T>
Matrix::Matrix<T> solveLaplace(
    const uint& size,
    const coord<T>& domain,
    const planeToScalar<T>& driver,
    const stencilGen<T>& makeStencil,
    const planeToScalar<T>& dirichlet
) {
  // Find limits of domain
  const auto a = std::get<0>(domain);
  const auto b = std::get<1>(domain);
  const auto h = (b - a) / (size - 1);

  // Internal limits of sampling
  const auto _a = a + h;
  const auto _b = a + b - h;
  const auto _intSize = size - 2;

  // Instantiate output grid (interior of mesh)
  Matrix::Matrix<T> rhs(
      _intSize,
      _intSize,
      [&](const uint& row, const uint& col) {
        // Determine location in domain
        const auto x = a + h * (col + 1);
        const auto y = a + h * (row + 1);

        // Initialize an accumulator
        T acc = 0;

        // Add on driving function
        acc += driver(x, y);

        // Determine stencil
        const auto samples = makeStencil({ x, y }, h);

        for (const auto& sample : samples) {
          const auto [ mult, coords ] = sample;
          const auto [ _x, _y ] = coords;

          // If outside of interior
          if (_x < _a || _x > _b || _y < _a || _y > _b) {
            // Subtract off the boundary (since we are truncating it)
            acc -= dirichlet(_x, _y);
          }
        }

        return acc;
      }
  );

  // Turn rhs into a column vector
  rhs = rhs.flatten();

  // Create the laplacian operator
  Matrix::Matrix<T> lap(_intSize * _intSize, _intSize * _intSize);

  // For each of the laplacian entries in lexigraphical order
  for (uint col = 0; col < _intSize; ++col) {
    for (uint row = 0; row < _intSize; ++row) {
      // Create a temporary matrix to embed stencil in
      Matrix::Matrix<T> stencilOp(_intSize, _intSize);

      // Determine stencil w.r.t. index, not domain
      const auto samples = makeStencil({ col, row }, 1);

      for (const auto& sample : samples) {
        const auto [ mult, coords ] = sample;
        const auto [ _col, _row ] = coords;

        // If inside of interior
        if (_col >= 0 && _col < _intSize && _row >= 0 && _row < _intSize) {
          // Embed point inside of matrix
          stencilOp.setVal(_col, _row, mult);
        }
      }

      // Flatten into a column vector
      stencilOp = stencilOp.flatten();

      // Add values to the correct row in laplace operator
      for (uint j = 0; j < _intSize * _intSize; ++j) {
        lap.setVal(row + _intSize * col, j, stencilOp.getVal(j, 0));
      }
    }
  }

  // Solve system for solution
  auto u = Matrix::Matrix<T>::solve(lap, rhs, Matrix::Solve::LU);

  // Reform matrix
  u = u.squareUp(_intSize, _intSize);

  return u;
}

int main() {
  // Define domain
  const coord<double> domain = { 0, 1 };

  // Define driving function
  const planeToScalar<double> driver =
    [](const double& x, const double& y) {
      return std::sin(x * y);
    };

  // G(x, y) along the boundaries
  const planeToScalar<double> boundary =
    [](const double& x, const double& y) -> double {
      (void) x;
      (void) y;

      return 0.0;
    };

  // Define size of mesh
  const uint size = 5;

  // Stencil generation
  const stencilGen<double> fivePoint = [](
      const coord<double>& center,
      const double& h
  ) {
    const auto [ x, y ] = center;
    const auto mult = 1.0 / (double) (h * h);

    return stencil<double>({
      { mult * -4,{ x + 0, y + 0 } },
      { mult * 1, { x + h, y + 0 } },
      { mult * 1, { x - h, y + 0 } },
      { mult * 1, { x + 0, y + h } },
      { mult * 1, { x + 0, y - h } }
    });
  };

  auto soln = solveLaplace<double>(size, domain, driver, fivePoint, boundary);

  std::cout << "Solution with 5-point stencil:" << std::endl;
  std::cout << soln << std::endl;

  // Output to file for funsies
  Image::ImageWriter::matrixHeatmap("soln.ppm", soln);

  return EXIT_SUCCESS;
}
