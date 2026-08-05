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
    return x >= number(0) ? number(1) : number(0);
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

  [[nodiscard]] RealType
  curvature(const Normal<dim, RealType> &n) const
  {
    if constexpr (dim == 2)
      {
        return {};
      }
    else if constexpr (dim == 3)
      {
        return curvature_3d(n);
      }
    else
      {
        UNREACHABLE("Anisotropies are only supported for 2D and 3D.");
      }
  }

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

    adouble a_x;
    a_x <<= n.x;
    adouble a_y = sqrt(RealType(1) - a_x * a_x);

    Normal<dim, adouble> a_n;
    a_n.x = a_x;
    a_n.y = a_y;

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

    Normal<dim, adouble> a_n;
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

private:
  Parameters<dim, RealType> param;
};

int
main(int argc, char *argv[])
{
  constexpr int dim = 3;
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

  RealType              min_K = std::numeric_limits<RealType>::max();
  Normal<dim, RealType> min_n;
  constexpr int         n_theta = 360;
  constexpr int         n_phi   = 180;
  for (int i = 0; i < n_theta; ++i)
    {
      const RealType theta = 2.0 * std::numbers::pi * i / n_theta;

      for (int j = 0; j <= n_phi; ++j)
        {
          const RealType        phi = std::numbers::pi * j / n_phi;
          Normal<dim, RealType> n(theta, phi);
          const auto            K = gamma.curvature(n);
          if (K < min_K)
            {
              min_K = K;
              min_n = n;
            }
        }
    }
  std::cout << "Minimum K = " << min_K << '\n'
            << "at (" << min_n.x << ", " << min_n.y << ", " << min_n.z << ")\n";

  return 0;
}
