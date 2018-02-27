---
math: true
permalink: /laplace
title: Solving Laplace Equation With FDM
layout: page
---

**Routine Name**: Laplace Equation With FDM

**Author**: Kyle Hovey

**Language**: C++

**Description/Purpose**:

This function solves the laplace equation over a square domain with Dirichlet boundary conditions that can be arbitrarily specified (homogeneous or not). In particular, we are interested in solving

# \\[ \nabla^2 u(x, y) = f(x, y) | (x, y) \in (a, b) \times (a, b) := D \\]

with Dirichlet boundary conditions

# \\[ g(x, y) | (x, y) \in \partial D \\]

**Input**:

Domain begin and end (`a` and `b`), size of mesh, function for computing stencil, and a function for the Dirichlet boundary conditions.

**Output**:

A 2D matrix that represents the solution along the interior of the domain.

**Usage/Example**:

First of all, I am using a few type conventions that will prove useful in making our definitions more terse and readable. They are as follows:

{% highlight C++ %}
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
{% endhighlight %}

The following code is inside of a rather large main function, so I'll break it up part by part to make the explanation more clear.

First, we generate the basic parameters that define the type of analysis we wish to perform. Here we define the domain, the function evaluated at the boundary, and the size of the mesh we wish to use.

{% highlight C++ %}
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
{% endhighlight %}

Then, we need to define what type of stencil we will be using for analysis. Here I define two stencils: the five and nine point equivalent discrete forms of the 2D laplace operator. Notice that the stencil points also include multipliers. To define a custom stencil, one must only create a new lambda function here and provide it to the solver.

{% highlight C++ %}
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

const stencilGen<double> ninePoint = [](
    const coord<double>& center,
    const double& h
) {
  const auto [ x, y ] = center;
  const auto mult = 1.0 / (double) (h * h);

  return stencil<double>({
    { mult * -8, { x - 0, y + 0 } },
    { mult * 1,  { x - 0, y + h } },
    { mult * 1,  { x - 0, y - h } },
    { mult * 1,  { x + h, y + 0 } },
    { mult * 1,  { x + h, y + h } },
    { mult * 1,  { x + h, y - h } },
    { mult * 1,  { x - h, y + 0 } },
    { mult * 1,  { x - h, y + h } },
    { mult * 1,  { x - h, y - h } }
  });
};
{% endhighlight %}

Then, to solve the equation, we just need to provide these parameters to the solver function.

{% highlight C++ %}
auto soln = solveLaplace<double>(size, domain, driver, fivePoint, boundary);

std::cout << "Solution with 5-point stencil:" << std::endl;
std::cout << soln << std::endl;

soln = solveLaplace<double>(size, domain, driver, ninePoint, boundary);

std::cout << "Solution with 9-point stencil:" << std::endl;
std::cout << soln << std::endl;
{% endhighlight %}

Output:

With a \\(25\\) point mesh (\\(9\\) interior points), we get the following result:

{% highlight C++ %}
Solution with 5-point stencil:
-0.0970489 -0.162868 -0.1537
-0.162868 -0.276049 -0.265528
-0.1537 -0.265528 -0.266089

Solution with 9-point stencil:
-0.0353471 -0.0593491 -0.0556861
-0.0593491 -0.101619 -0.0981173
-0.0556861 -0.0981173 -0.103895
{% endhighlight %}

This is much more interesting when we increase the mesh size to a \\(200 \times 200 \\) grid. Since it would be nearly incomprehensible to look at a solution that large with just textual output, I wrote a method to convert a matrix into a heatmap where the color range is restricted linearly between the min and max values of the matrix. Here is the aformentioned higher-res solution:

![Laplace Solution Image](/MATH_5620/images/laplace.jpg)

**Implementation/Code:**

To implement this solver in a general way, I use lambda functions as parameters for both the Dirichlet boundary conditions and the stencil desired to be used. Both the right-hand-side and discrete Laplace operator are generated on the fly using these parameters. I make use of an additional method I added to my matrix class wherein I may flatten a matrix into a column vector according to lexigraphical ordering, or reshape it into an \\(m\\) by \\(n\\) matrix if it represents a column vector. This way I can mutate a matrix in a way that makes sense, then flatten it into the ordering required for our linear system. The comments in my code should help explain further:

{% highlight c++ %}
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
  auto u = Matrix::Matrix<T>::solve(lap, rhs, Matrix::Solve::Jacobi);

  // Reform matrix
  u = u.squareUp(_intSize, _intSize);

  return u;
}
{% endhighlight %}
