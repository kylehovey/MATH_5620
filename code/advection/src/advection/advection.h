#ifndef ADVECTION_H
#define ADVECTION_H

#include <iostream>
#include "../../../matrix/src/matrix/matrix.h"
#include "../../../image/src/image/image.h"

namespace Advection {
  using uint = unsigned int;

  template <typename T>
  using interval = std::tuple<T, T>;

  template <typename T>
  using monad = std::function<T(const T&)>;

  /**
   * Test the upwinding method for the PDE:
   *  u_t + a u_x = 0; boundary conditions: 0
   * @param size Size of mesh
   * @param spaceDomain Interval in space to evaluate over
   * @param timeDomain Interval in time to evaluate over
   * @param spaceStep Space step quantity
   * @param timeStep Time step quantity
   * @param eta Initial conditions along t = t_0
   * @param a Value seen in advection equation
   */
  template<typename T>
  void testUpwinding(
      const uint size,
      const interval<T> spaceDomain,
      const interval<T> timeDomain,
      const T& spaceStep,
      const T& timeStep,
      const monad<T>& eta,
      const T& a
  ) {
    // Solution vector
    Matrix::Matrix<T> U(size - 2, 1, [=](const uint& a, const uint& b) -> T {
        (void) b;

        const auto A = std::get<0>(spaceDomain);
        const auto x = A + (a + 1) * spaceStep;

        return eta(x);
    });

    // Time-step operator
    const auto coeffs = (std::vector<T>) {
      a * spaceStep / timeStep,
      1 - a * spaceStep / timeStep
    };
    const auto A = Matrix::Matrix<T>::genNDiag(
        size - 2, 
        coeffs, 
        -1
    );

    // Evolve system in time
    std::vector<Matrix::Matrix<T>> solns = { U };
    const auto [ t_0, t_f ] = timeDomain;
    for (auto t = t_0; t < t_f; t += timeStep) {
      U = A * U;
      solns.push_back(U);
    }

    const Matrix::Matrix<T> wholePicture(
        size - 2,
        solns.size(),
        [=](const uint& a, const uint& b) -> double {
          return solns[b].getVal(a, 0);
        }
    );

    // Output
    Image::ImageWriter::matrixHeatmap("./upWinding.ppm", wholePicture, 1000);
  }

  /**
   * Test the Lax-Wendroff method for the PDE:
   *  u_t + a u_x = 0; boundary conditions: 0
   * @param size Size of mesh
   * @param spaceDomain Interval in space to evaluate over
   * @param timeDomain Interval in time to evaluate over
   * @param spaceStep Space step quantity
   * @param timeStep Time step quantity
   * @param eta Initial conditions along t = t_0
   * @param a Value seen in advection equation
   */
  template<typename T>
  void testLaxWendroff(
      const uint size,
      const interval<T> spaceDomain,
      const interval<T> timeDomain,
      const T& spaceStep,
      const T& timeStep,
      const monad<T>& eta,
      const T& a
  ) {
    // Solution vector
    Matrix::Matrix<T> U(size - 2, 1, [=](const uint& a, const uint& b) -> T {
        (void) b;

        const auto A = std::get<0>(spaceDomain);
        const auto x = A + (a + 1) * spaceStep;

        return eta(x);
    });

    // Time-step operator
    const auto alpha = a * spaceStep / (2 * timeStep);
    const auto beta = a * a * spaceStep * spaceStep / (2 * timeStep * timeStep);
    const auto coeffs = (std::vector<T>) {
      alpha + beta,   // U_{j-1}
      1 - 2 * beta,   // U_{j}
      -alpha + beta,  // U_{j+1}
    };
    const auto A = Matrix::Matrix<T>::genNDiag(
        size - 2, 
        coeffs, 
        -1
    );

    // Evolve system in time
    std::vector<Matrix::Matrix<T>> solns = { U };
    const auto [ t_0, t_f ] = timeDomain;
    for (auto t = t_0; t < t_f; t += timeStep) {
      U = A * U;
      solns.push_back(U);
    }

    const Matrix::Matrix<T> wholePicture(
        size - 2,
        solns.size(),
        [=](const uint& a, const uint& b) -> double {
          return solns[b].getVal(a, 0);
        }
    );

    // Output
    Image::ImageWriter::matrixHeatmap("./laxWendroff.ppm", wholePicture);
  }

  /**
   * Test the Beam Warming method for the PDE:
   *  u_t + a u_x = 0; boundary conditions: 0
   * @param size Size of mesh
   * @param spaceDomain Interval in space to evaluate over
   * @param timeDomain Interval in time to evaluate over
   * @param spaceStep Space step quantity
   * @param timeStep Time step quantity
   * @param eta Initial conditions along t = t_0
   * @param a Value seen in advection equation
   */
  template<typename T>
  void testBeamWarming(
      const uint size,
      const interval<T> spaceDomain,
      const interval<T> timeDomain,
      const T& spaceStep,
      const T& timeStep,
      const monad<T>& eta,
      const T& a
  ) {
    // Solution vector
    Matrix::Matrix<T> U(size - 2, 1, [=](const uint& a, const uint& b) -> T {
        (void) b;

        const auto A = std::get<0>(spaceDomain);
        const auto x = A + (a + 1) * spaceStep;

        return eta(x);
    });

    // Time-step operator
    const auto alpha = a * spaceStep / (2 * timeStep);
    const auto beta = a * a * spaceStep * spaceStep / (2 * timeStep * timeStep);
    const auto coeffs = (std::vector<T>) {
      beta - alpha,
      4 * alpha - 2 * beta,
      1 - 3 * alpha + beta
    };
    const auto A = Matrix::Matrix<T>::genNDiag(
        size - 2, 
        coeffs, 
        -2
    );

    // Evolve system in time
    std::vector<Matrix::Matrix<T>> solns = { U };
    const auto [ t_0, t_f ] = timeDomain;
    for (auto t = t_0; t < t_f; t += timeStep) {
      U = A * U;
      solns.push_back(U);
    }

    const Matrix::Matrix<T> wholePicture(
        size - 2,
        solns.size(),
        [=](const uint& a, const uint& b) -> double {
          return solns[b].getVal(a, 0);
        }
    );

    // Output
    Image::ImageWriter::matrixHeatmap("./beamWarming.ppm", wholePicture);
  }
}

#endif
