#include <iostream>
#include "./advection/advection.h"

int main() {
  const auto spaceDomain = std::make_tuple(0.0, 1.0);
  const auto timeDomain = std::make_tuple(0.0, 1.0);
  const auto dt = 1e-3;
  const auto dx = 1e-3;
  const auto [ A, B ] = spaceDomain;
  const auto size = std::round((B - A) / dx);
  const auto a = 0.7;
  const auto eta = [](const double& x) -> double {
    return (x >= 0.3 && x <= 0.6) ? 100 : 0;
  };

  Advection::testUpwinding<double>(
      size,
      spaceDomain,
      timeDomain,
      dx,
      dt,
      eta,
      a
  );

  Advection::testLaxWendroff<double>(
      size,
      spaceDomain,
      timeDomain,
      dx,
      dt,
      eta,
      a
  );

  Advection::testBeamWarming<double>(
      size,
      spaceDomain,
      timeDomain,
      dx,
      dt,
      eta,
      a
  );

  return EXIT_SUCCESS;
}
