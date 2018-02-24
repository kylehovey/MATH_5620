#include "./image/image.h"
#include <iostream>

int main() {
  const auto A = Matrix::Matrix<double>::hilbert(100);

  Image::ImageWriter::matrixHeatmap("out.ppm", A);

  std::cout << A << std::endl;

  return EXIT_SUCCESS;
}
