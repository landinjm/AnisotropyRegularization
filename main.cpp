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
    return surface_energy(param, n);
  }

  /**
   * Heaviside step function
   */
  [[nodiscard]] constexpr RealType
  theta(RealType x) const
  {
    return x >= RealType {0} ? RealType {1} : RealType {0};
  }

  /**
   * Surface energy
   */
  [[nodiscard]] RealType
  surface_energy(const Parameters<dim, RealType> &param,
                 const Normal<dim, RealType>     &n) const
  {
    auto sum_term = RealType {0};
    for (int i = 0; i < param.N; i++)
      {
        const auto x = n.dot(param.m[i]);
        sum_term += param.alpha[i] * std::pow(x, param.omega[i]) * theta(x);
      }
    return param.gamma_0 * (RealType {1} - sum_term);
  }

private:
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

  return 0;
}
