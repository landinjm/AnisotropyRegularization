#include <cmath>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <libassert/assert.hpp>
#include <numbers>
#include <string>
#include <vector>

template <int dim, typename RealType>
struct Normal;

template <typename RealType>
struct Normal<2, RealType>
{
  RealType x;
  RealType y;

  Normal()
    : x(RealType {1})
    , y(RealType {0})
  {}

  explicit Normal(RealType theta)
    : x(std::cos(theta))
    , y(std::sin(theta))
  {}

  Normal(RealType x_, RealType y_)
  {
    RealType d = std::sqrt(x_ * x_ + y_ * y_);
    DEBUG_ASSERT(d > RealType {0}, "Vector input length must be nonzero");

    x = x_ / d;
    y = y_ / d;
  }

  constexpr RealType
  dot(const Normal &other) const
  {
    return x * other.x + y * other.y;
  }
};

template <typename RealType>
struct Normal<3, RealType>
{
  RealType x;
  RealType y;
  RealType z;

  Normal()
    : x(RealType {1})
    , y(RealType {0})
    , z(RealType {0})
  {}

  Normal(RealType theta, RealType psi)
    : x(std::cos(psi) * std::cos(theta))
    , y(std::cos(psi) * std::sin(theta))
    , z(std::sin(psi))
  {}

  Normal(RealType x_, RealType y_, RealType z_)
  {
    RealType d = std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    DEBUG_ASSERT(d > RealType {0}, "Vector input length must be nonzero");

    x = x_ / d;
    y = y_ / d;
    z = z_ / d;
  }

  constexpr RealType
  dot(const Normal &other) const
  {
    return x * other.x + y * other.y + z * other.z;
  }
};

template <int dim, typename RealType>
struct Parameters
{
  /**
   * Baseline surface energy
   */
  RealType gamma_0;

  /**
   * Total number of energy minima
   */
  unsigned int N;

  /**
   * Unit vectors that give the energy minma
   */
  std::vector<Normal<dim, RealType>> m;

  /**
   * Energy depth of the minima
   */
  std::vector<RealType> alpha;

  /**
   * Energy width of the minima
   */
  std::vector<RealType> omega;

  /**
   * Resize the depth and width of minima to the total number of minima
   */
  void
  resize()
  {
    m.resize(N);
    alpha.resize(N);
    omega.resize(N);
  }

  /**
   * Validate the parameters
   */
  void
  validate() const
  {
    DEBUG_ASSERT(N >= 2, "There must be at least two energy minima");
    DEBUG_ASSERT(
      m.size() == N && alpha.size() == N && omega.size() == N,
      "The direction, depth, and width values of the energy minima do not match the "
      "number of energy minima. Make sure to call the `resize()` function.",
      N,
      m.size(),
      alpha.size(),
      omega.size());
  }
};

/**
 * Heaviside step function
 */
template <typename RealType>
inline constexpr RealType
theta(RealType x)
{
  return x >= RealType {0} ? RealType {1} : RealType {0};
}

/**
 * Surface energy
 */
template <int dim, typename RealType>
inline RealType
surface_energy(const Parameters<dim, RealType> &param, Normal<dim, RealType> n)
{
  auto sum_term = RealType {0};
  for (int i = 0; i < param.N; i++)
    {
      const auto x = n.dot(param.m[i]);
      sum_term += param.alpha[i] * std::pow(x, param.omega[i]) * theta(x);
    }
  return param.gamma_0 * (RealType {1} - sum_term);
}

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

  std::string keyword;

  prm >> keyword;
  DEBUG_ASSERT(keyword == "gamma0", "Expected keyword gamma0", keyword);
  prm >> param.gamma_0;

  prm >> keyword;
  DEBUG_ASSERT(keyword == "N", "Expected keyword N", keyword);
  prm >> param.N;

  param.resize();

  prm >> keyword;
  DEBUG_ASSERT(keyword == "m", "Expected keyword m", keyword);
  for (int i = 0; i < param.N; ++i)
    {
      RealType x {};
      RealType y {};
      RealType z {};
      prm >> x >> y >> z;
      param.m[i] = Normal<dim, RealType>(x, y, z);
    }

  prm >> keyword;
  DEBUG_ASSERT(keyword == "alpha", "Expected keyword alpha", keyword);
  for (int i = 0; i < param.N; ++i)
    {
      prm >> param.alpha[i];
    }

  prm >> keyword;
  DEBUG_ASSERT(keyword == "omega", "Expected keyword omega", keyword);
  for (int i = 0; i < param.N; ++i)
    {
      prm >> param.omega[i];
    }

  param.validate();

  const int Ntheta = 181;
  const int Npsi   = 91;

  std::ofstream vtk("surface.vtk");
  DEBUG_ASSERT(vtk, "Could not open output file.");

  vtk << "# vtk DataFile Version 3.0\n";
  vtk << "Surface energy\n";
  vtk << "ASCII\n";
  vtk << "DATASET STRUCTURED_GRID\n";
  vtk << "DIMENSIONS " << Ntheta << " " << Npsi << " 1\n";
  vtk << "POINTS " << Ntheta * Npsi << " double\n";

  std::vector<RealType> gamma_values;
  gamma_values.reserve(Ntheta * Npsi);

  for (int j = 0; j < Npsi; ++j)
    {
      RealType psi = -M_PI / 2 + M_PI * j / (Npsi - 1);

      for (int i = 0; i < Ntheta; ++i)
        {
          RealType theta = 2 * M_PI * i / (Ntheta - 1);

          Normal<dim, RealType> n(theta, psi);

          RealType gamma = surface_energy(param, n);

          gamma_values.push_back(gamma);

          vtk << gamma * n.x << ' ' << gamma * n.y << ' ' << gamma * n.z << '\n';
        }
    }

  vtk << "\nPOINT_DATA " << Ntheta * Npsi << '\n';
  vtk << "SCALARS gamma double\n";
  vtk << "LOOKUP_TABLE default\n";

  for (RealType g : gamma_values)
    vtk << g << '\n';

  return 0;
}
