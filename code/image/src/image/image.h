#ifndef IMAGE_H
#define IMAGE H

#include "../../../matrix/src/matrix/matrix.h"
#include <fstream>
#include <cmath>

namespace Image {
  struct Color {
    Color(
        const short int& R,
        const short int& G,
        const short int& B
        ) : R(R), G(G), B(B) { };

    short int R, G, B;
  };

  const Color Red   (255, 0, 0);
  const Color Green (0, 255, 0);
  const Color Blue  (0, 0, 255);

  std::function<Color(double)> lerp(const Color& a, const Color& b) {
    return [=](const double& t) {
      const auto s = 1 - t;

      return Color(
          s * a.R + t * b.R,
          s * a.G + t * b.G,
          s * a.B + t * b.B
      );
    };
  }

  class ImageWriter {
    public:
      /**
       * Write a matrix to a ppm file
       * @param path File path to write to
       * @param grid Matrix of values to write out
       * @param width Width of output image in pixels
       */
      template <typename T>
      static void matrixHeatmap(
          const std::string& path,
          const Matrix::Matrix<T>& grid,
          const uint& width = 1000
      );
    private:
      ImageWriter(void) {}
  };

  template <typename T>
  void ImageWriter::matrixHeatmap(
      const std::string& path,
      const Matrix::Matrix<T>& grid,
      const uint& width
  ) {
    std::ofstream outFile;
    outFile.open(path);

    // RGB format
    outFile << "P3" << std::endl;

    // Dimension of output
    const auto [ m, n ] = grid.getSize();
    const double aspect = (double) n / m;
    const uint height = std::round((double) width / aspect);
    const double ppcX = (double) width / n;
    const double ppcY = (double) width / m;

    outFile << width << " " << height << std::endl;

    // Color resolution
    outFile << 255 << std::endl;

    // Determine range of values
    const auto min = grid.getMin();
    const auto range = grid.getMax() - min;

    // Color range
    const auto wheel = Image::lerp(Image::Blue, Image::Red);

    // For each pixel
    for (uint y = 0; y < height; ++y) {
      for (uint x = 0; x < width; ++x) {
        // Determine cell of matrix this pixel is iside of
        const auto row = std::floor((double) y / ppcY);
        const auto col = std::floor((double) x / ppcX);

        // Get value in cell
        const auto val = grid.getVal(row, col);

        // Get color for cell
        const auto color = wheel((double) (val - min) / range);

        // Write value out
        outFile << color.R << " ";
        outFile << color.G << " ";
        outFile << color.B << "   ";
      }

      outFile << std::endl;
    }

    outFile.close();
  }
};

#endif
