#include "src/normal.h"
#include "src/parameters.h"
#include "src/vtk.h"

#include <adolc/adolc.h>
#include <cmath>
#include <fstream>
#include <libassert/assert.hpp>

template <int dim, typename RealType>
class SurfaceEnergy
{
public:
  explicit SurfaceEnergy(const Parameters<dim, RealType> &_param)
    : param(_param)
  {}

  RealType
  operator()(const Normal<dim, RealType> &n) const
  {
    return surface_energy<RealType>(n);
  }

  /**
   * Heaviside step function
   *
   * @note Allows us to control the positive and negative directions separately
   */
  template <typename number>
  [[nodiscard]] constexpr number
  theta(number x) const
  {
    using std::tanh;

    constexpr double eps = 1.0e-6;
    return number(0.5) * (number(1) + tanh(x / number(eps)));
  }

  /**
   * Surface energy
   */
  template <typename number>
  [[nodiscard]] number
  surface_energy(const Normal<dim, number> &n) const
  {
    using std::pow;

    auto sum_term = number(0);
    for (int i = 0; i < param.N; i++)
      {
        const auto x = n.dot(param.m[i]);
        sum_term += param.alpha[i] * pow(x, param.omega[i]) * theta(x);
      }
    return param.gamma_0 * (number(1) - sum_term);
  }

  [[nodiscard]] Normal<dim, RealType>
  wulff_point(const Normal<dim, RealType> &n) const
  {
    if constexpr (dim == 2)
      {
        return wulff_point_2d(n);
      }
    else if constexpr (dim == 3)
      {
        return wulff_point_3d(n);
      }
    else
      {
        static_assert(dim == 2 || dim == 3,
                      "Anisotropies are only supported for 2D and 3D");
      }
  }

  [[nodiscard]] Normal<dim, RealType>
  wulff_point_2d(const Normal<dim, RealType> &n) const
  {
    // NOLINTBEGIN

    // Convert normal into array that adolc can use.
    RealType theta  = std::atan2(n.y, n.x);
    double   x_0[1] = {theta};

    // Tape the evaluation of the surface energy function
    const short tag = 22;
    trace_on(tag);

    adouble a_theta;

    a_theta <<= theta;

    Normal<2, adouble> a_n(a_theta);

    adouble a_gamma = surface_energy(a_n);

    double gamma;
    a_gamma >>= gamma;

    trace_off();

    // Get the gradient of gamma
    // 1 independent variables (theta)
    double grad[1];
    gradient(tag, 1, x_0, grad);

    // Convert to RealType
    const RealType gamma_theta = grad[0];

    // Get the tangent
    Normal<2, RealType> t;
    t.x = -std::sin(theta);
    t.y = std::cos(theta);

    // Compute the wulff tangent
    Normal<2, RealType> x;
    x.x = gamma * n.x + gamma_theta * t.x;
    x.y = gamma * n.y + gamma_theta * t.y;

    return x;

    // NOLINTEND
  }

  [[nodiscard]] Normal<dim, RealType>
  wulff_point_3d(const Normal<dim, RealType> &n) const
  {
    // NOLINTBEGIN

    // Convert normal into array that adolc can use.
    RealType theta  = std::atan2(n.y, n.x);
    RealType psi    = std::asin(n.z);
    double   x_0[2] = {theta, psi};

    // Tape the evaluation of the surface energy function
    const short tag = 33;
    trace_on(tag);

    adouble a_theta;
    adouble a_psi;

    a_theta <<= theta;
    a_psi <<= psi;

    Normal<3, adouble> a_n(a_theta, a_psi);

    adouble a_gamma = surface_energy(a_n);

    double gamma;
    a_gamma >>= gamma;

    trace_off();

    // Get the gradient of gamma
    // 2 independent variables (theta, psi)
    double grad[2];
    gradient(tag, 2, x_0, grad);

    // Convert to RealType
    const RealType gamma_theta = grad[0];
    const RealType gamma_psi   = grad[1];

    // Get the tangent
    Normal<3, RealType> e_theta;
    e_theta.x = -std::sin(theta);
    e_theta.y = std::cos(theta);
    e_theta.z = 0.0;

    Normal<3, RealType> e_psi;
    e_psi.x = -std::cos(theta) * std::sin(psi);
    e_psi.y = -std::sin(theta) * std::sin(psi);
    e_psi.z = std::cos(psi);

    // Compute the wulff tangent
    Normal<3, RealType> x;

    const RealType c = std::cos(psi);

    x.x = gamma * n.x + gamma_theta * e_theta.x / c + gamma_psi * e_psi.x;
    x.y = gamma * n.y + gamma_theta * e_theta.y / c + gamma_psi * e_psi.y;
    x.z = gamma * n.z + gamma_psi * e_psi.z;

    return x;

    // NOLINTEND
  }

  /**
   * Curvature
   */
  [[nodiscard]] RealType
  curvature(const Normal<dim, RealType> &n) const
  {
    if constexpr (dim == 2)
      {
        return curvature_2d(n);
      }
    else if constexpr (dim == 3)
      {
        return curvature_3d(n);
      }
    else
      {
        static_assert(dim == 2 || dim == 3,
                      "Anisotropies are only supported for 2D and 3D");
      }
  }

  bool
  requires_regularization(int n_theta = 181, int n_psi = 91) const
  {
    RealType              min_curvature = std::numeric_limits<RealType>::max();
    Normal<dim, RealType> min_n;

    if constexpr (dim == 2)
      {
        for (int i = 0; i < n_theta; ++i)
          {
            RealType              theta = 2.0 * std::numbers::pi * i / n_theta;
            Normal<dim, RealType> n(theta);
            const auto            K = curvature_2d(n);
            if (K < min_curvature)
              {
                min_curvature = K;
                min_n         = n;
              }
          }
        std::cout << "Minimum curvature = " << min_curvature << '\n'
                  << "at (" << min_n.x << ", " << min_n.y << ")\n";
      }
    else if constexpr (dim == 3)
      {
        for (int j = 0; j < n_psi; ++j)
          {
            RealType psi = -std::numbers::pi / 2.0 + std::numbers::pi * j / n_psi;
            for (int i = 0; i < n_theta; ++i)
              {
                RealType              theta = 2.0 * std::numbers::pi * i / n_theta;
                Normal<dim, RealType> n(theta, psi);
                const auto            K = curvature_3d(n);
                if (K < min_curvature)
                  {
                    min_curvature = K;
                    min_n         = n;
                  }
              }
          }
        std::cout << "Minimum curvature = " << min_curvature << '\n'
                  << "at (" << min_n.x << ", " << min_n.y << ", " << min_n.z << ")\n";
      }
    else
      {
        static_assert(dim == 2 || dim == 3,
                      "Anisotropies are only supported for 2D and 3D");
      }

    return min_curvature < 0.0;
  }

private:
  [[nodiscard]] RealType
  curvature_2d(const Normal<dim, RealType> &n) const
  {
    // NOLINTBEGIN

    // Convert normal into array that adolc can use. For simplicity, we just differentiate
    // with respect to theta
    RealType theta  = std::atan2(n.y, n.x);
    double   x_0[1] = {theta};

    // Tape the evaluation of the surface energy function
    const short tag = 2;
    trace_on(tag);

    adouble a_theta;

    a_theta <<= theta;

    Normal<2, adouble> a_n(a_theta);

    adouble a_gamma = surface_energy(a_n);

    double gamma;
    a_gamma >>= gamma;

    trace_off();

    // Get the hessian of gamma
    // 1 independent variables (theta)
    double  hess[1][1];
    double *hess_ptr[1] = {hess[0]};
    hessian(tag, 1, x_0, hess_ptr);

    // Convert to RealType
    // NOTE: The hessian only returns the lower triangle
    const RealType gamma_theta_theta = hess[0][0];

    // Evaluate curvature
    return RealType(1) / (gamma + gamma_theta_theta);
    // NOLINTEND
  }

  [[nodiscard]] RealType
  curvature_3d(const Normal<dim, RealType> &n) const
  {
    // NOLINTBEGIN

    // Convert normal into array that adolc can use. n_z is constraint by the other two so
    // we don't need to differentiate with respect to it.
    double x_0[2] = {n.x, n.y};

    // Tape the evaluation of the surface energy function
    const short tag = 3;
    trace_on(tag);

    adouble a_x;
    adouble a_y;
    a_x <<= n.x;
    a_y <<= n.y;
    adouble a_z = sqrt(RealType(1) - a_x * a_x - a_y * a_y);

    Normal<3, adouble> a_n;
    a_n.x = a_x;
    a_n.y = a_y;
    a_n.z = a_z;

    adouble a_gamma = surface_energy(a_n);

    double gamma;
    a_gamma >>= gamma;

    trace_off();

    // Get gradient of gamma
    // 2 independent variables (n_x, n_y)
    double grad[2];
    gradient(tag, 2, x_0, grad);

    // Get the hessian of gamma
    // 2 independent variables (n_x, n_y)
    double  hess[2][2];
    double *hess_ptr[2] = {hess[0], hess[1]};
    hessian(tag, 2, x_0, hess_ptr);

    // Convert to RealType
    // NOTE: The hessian only returns the lower triangle
    const RealType gamma_x = grad[0];
    const RealType gamma_y = grad[1];

    const RealType gamma_xx = hess[0][0];
    const RealType gamma_xy = hess[1][0];
    const RealType gamma_yy = hess[1][1];

    // Evaluate the Gaussian curvature
    const RealType v = gamma - n.x * gamma_x - n.y * gamma_y;
    const RealType u = gamma_xx * (RealType(1) - n.x * n.x) +
                       gamma_yy * (RealType(1) - n.y * n.y) -
                       RealType(2) * n.x * n.y * gamma_xy;

    const RealType g_x = v * n.x + gamma_x;
    const RealType g_y = v * n.y + gamma_y;
    const RealType g_z = v * n.z;

    const RealType xi = std::sqrt(g_x * g_x + g_y * g_y + g_z * g_z);
    const RealType r  = gamma / xi;
    const RealType h  = gamma_xx * gamma_yy - gamma_xy * gamma_xy;

    return (r * r * r * r) * (v * v + v * u + n.z * n.z * h);
    // NOLINTEND
  }

  Parameters<dim, RealType> param;
};

int
main(int argc, char *argv[])
{
  constexpr int dim = 2;
  using RealType    = double;

  // Try to open the parameter file
  std::string   filename = (argc > 1) ? argv[1] : "parameter.prm";
  std::ifstream prm(filename);
  DEBUG_ASSERT(prm, "Could not open parameter file.", filename);

  // Read the parameter inputs
  Parameters<dim, RealType> param;
  param.read(prm);
  param.validate();

  // Output the surface energy surface to vtk
  SurfaceEnergy<dim, RealType> gamma(param);
  write_gamma_surface_vtk<dim, RealType>("surface.vtk", gamma);

  // Check curvature
  if (gamma.requires_regularization(360, 180))
    {
      std::cout
        << "Warning: These values for the interfacial energy require regularization!"
        << std::endl;
    }

  write_wulff_surface_vtk<dim, RealType>("wulff.vtk", gamma);

  return 0;
}
