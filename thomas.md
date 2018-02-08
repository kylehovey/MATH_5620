---
math: true
permalink: /thomas
title: Thomas Algorithm, Elliptic ODE, and Error Vector Norm
layout: page
---

**Routine Name**: Thomas Algorithm, Elliptic ODE, and Error Vector Norm

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This algorithm uses a finite difference approximation of a second order derivative to solve the equation:

# \\[ u'' = f \\]

This equation has solution:

# \\[ u(x) = c_1 x + c_2 + \frac{\sin(\pi x)}{\pi^2} \\]

For the interval \\( (0, 1) \\), initial conditions \\( u(0) = 2.5, u(1) = 5 \\), and driving function \\( f(x) = \sin(\pi x) \\), the exact solution becomes:

# \\[ u(x) = 2.5x + 2.5 + \frac{\sin(\pi x)}{\pi^2} \\]

**Input**:

The interval \\( \(a, b\) \\), \\( u(a), u(b) \\), the driving function \\( f \\), and the size of the mesh.

**Output**:

A column vector representing the approximation of the solution \\( u \\) at the mesh points.

**Usage/Example**:

{% highlight C++ %}
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
  Mtx uExact(meshSize - 1, 1, [&](const uint& i, const uint& j) {
      (void) j;
      return exact(a + (i + 1) * h);
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
{% endhighlight %}

Output:

{% highlight C++ %}
Solved solution
2.71843
2.93995
3.16735
3.40284
3.64784
3.90284
4.16735
4.43995
4.71843

Exact solution
2.71869
2.94044
3.16803
3.40364
3.64868
3.90364
4.16803
4.44044
4.71869

Error vector
-0.00025879
-0.000492248
-0.000677521
-0.000796474
-0.000837462
-0.000796474
-0.000677521
-0.000492248
-0.00025879

1-norm of error vector
0.00528753
2-norm of error vector
0.00187262
infinity-norm of error vector
-0.00025879
{% endhighlight %}

**Implementation/Code:**

From the first problem on HW 2, I created a method that solves for finite difference coefficients of arbitrary order. To solve this code, I created a function that puts those coefficients into a matrix such that the coefficients are staggered across the diagonal of the matrix. For instance, we need a second order derivative here which yields a tri-diagonal matrix with \\(1, -2, 1\\) as the coefficients across the diagonal.

To handle boundary conditions, I subtract them off at the boundaries (as seen in the textbook). Then, to see the error of the method, I compute the exact solution over the mesh chosen and take the 2-norm of the difference between the exact solution and the approximated solution. I do this using the Thomas Algorithm, which implements a linear scan across the diagonal. Since the matrix is of dimension \\(n \times n\\), then the diagonal is of size \\(n\\) and the Thomas Algorithm solves this system in:

# \\[ O(n) \\]

Following is the code required to complete these goals:

{% highlight c++ %}
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
  const Matrix::Matrix<T> F(n - 1, 1, [&](const uint& i, const uint& j) -> T {
      (void) j;
      return std::pow(h, 2) * f(a + (i + 1) * h);
  });

  // Use initial conditions
  Matrix::Matrix<T> init(n - 1, 1);
  init.setVal(0, 0, ua);
  init.setVal(n - 2, 0, ub);

  // Generate differential operator for second derivative
  const auto D = Matrix::Matrix<T>::genFDMatrix(n - 1, 2);

  // Solve for solution vector
  const auto soln = Matrix::Matrix<T>::solve(
      D,
      F - init,
      Matrix::Solve::Thompson
   );

  return soln;
}
{% endhighlight %}

Where the thomas algorithm is implemented as follows (this is just the Thomas part of the general `solve` function):

{% highlight c++ %}
template <typename T>
Matrix<T> Matrix<T>::solve(
    const Matrix<T>& A,
    const Matrix<T>& b
) {
  // Only work for square matrices
  if (!A.isSquare()) {
    throw std::domain_error("Input matrix must be square.");
  }

  // Get size of matrix
  const auto m = std::get<0>(A.getSize());

  if (!A.isNDiagonal(3)) {
    throw std::domain_error("Thomas method needs tri-diagonal matrix.");
  }

  // Throw diagonals into vectors (store in a vector of diags)
  std::vector<std::vector<T>> diags(3, std::vector<T>());

  diags[0].push_back(0);

  for (uint row = 0; row < m; ++row) {
    for (uint nDiag = 0; nDiag < 3; ++nDiag) {
      const auto pos = row - 1 + nDiag;

      if (pos < m) {
        diags[nDiag].push_back(A.getVal(row, (uint) pos));
      }
    }
  }

  diags[2].push_back(0);

  // Syntactic nicety
  auto _a = diags[0];
  auto _b = diags[1];
  auto _c = diags[2];

  // Copy solution vector
  auto _d = std::vector<T>();
  for (uint i = 0; i < m; ++i) {
    _d.push_back(b.getVal(i, 0));
  }

  // Do elimination via algebra
  _c[0] = _c[0] / _b[0];
  _d[0] = _d[0] / _b[0];
  for (uint i = 1; i < m; ++i) {
    _c[i] = _c[i] / (_b[i] - _a[i] * _c[i - 1]);
    _d[i] = (_d[i] - _a[i] * _d[i - 1]) / (_b[i] - _a[i] * _c[i - 1]);
  }

  // Create the solution
  Matrix<T> x(m, 1);
  x.setVal(m - 1, 0, _d[m - 1]);
  for (int i = m - 2; i >= 0; --i) {
    const uint _i = i;

    x.setVal(_i, 0, _d[_i] - _c[_i] * x.getVal(_i + 1, 0));
  }

  return x;
}
{% endhighlight %}

The vector norm is calculated using a general function that will calculate any \\( p \\)-norm for a vector of arbitrary size. If \\( p \\) is the maximum value for its type, then the infinity norm is given instead.

{% highlight c++ %}
template <typename T>
T Matrix<T>::vNorm(const Matrix<T>& v, const uint& n) {
  const auto [ M, N ] = v.getSize();
  if (M != 1 && N != 1 && M != N) {
    throw std::domain_error("Need a row or column vector for vector norm.");
  }

  const auto isRow = (M == 1);
  const auto size = isRow ? N : M;

  T sum = 0;
  T max = std::numeric_limits<T>::min();

  for (uint i = 0; i < size; ++i) {
    const auto val = isRow ? v.getVal(0, i) : v.getVal(i, 0);

    max = (val > max) ? val : max;
    sum += pow(std::abs(val), n);
  }

  if (n == std::numeric_limits<uint>::max()) {
    // Infinity norm
    return max;
  } else {
    return pow(sum, 1.0 / n);
  }
}
{% endhighlight %}
