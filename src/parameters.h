#pragma once

#include "normal.h"

#include <fstream>
#include <libassert/assert.hpp>
#include <sstream>
#include <string>
#include <vector>

template <int dim, typename RealType>
struct Parameters
{
  static_assert(dim == 2 || dim == 3, "Anisotropies are only supported for 2D and 3D");

  /**
   * Baseline surface energy
   */
  RealType gamma_0 {};

  /**
   * Total number of energy minima
   */
  unsigned int N {};

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
   * Read parameters from input stream
   */
  void
  read(std::ifstream &input_stream)
  {
    std::string keyword;

    input_stream >> keyword;
    DEBUG_ASSERT(keyword == "gamma0", "Expected keyword gamma0", keyword);
    input_stream >> gamma_0;

    input_stream >> keyword;
    DEBUG_ASSERT(keyword == "N", "Expected keyword N", keyword);
    input_stream >> N;

    resize();

    input_stream >> keyword;
    DEBUG_ASSERT(keyword == "m", "Expected keyword m", keyword);
    for (int i = 0; i < N; ++i)
      {
        RealType x {};
        RealType y {};
        RealType z {};
        if constexpr (dim == 2)
          {
            input_stream >> x >> y;
            m[i] = Normal<2, RealType>(x, y);
          }
        else if constexpr (dim == 3)
          {
            input_stream >> x >> y >> z;
            m[i] = Normal<3, RealType>(x, y, z);
          }
        else
          {
            UNREACHABLE("Anisotropies are only supported for 2D and 3D.");
          }
      }

    // For this next bit we get the depths and widths. These can come in two flavors: a
    // single constant value and N values.
    input_stream >> keyword;
    DEBUG_ASSERT(keyword == "alpha", "Expected keyword alpha", keyword);

    std::string line;
    std::getline(input_stream >> std::ws, line);
    std::istringstream    iss(line);
    std::vector<RealType> values;
    RealType              value {};
    while (iss >> value)
      {
        values.push_back(value);
      }
    if (values.size() == 1)
      {
        std::fill(alpha.begin(), alpha.end(), values[0]);
      }
    else
      {
        DEBUG_ASSERT(values.size() == N, "Expected 1 or N alpha values", values.size());
        alpha = values;
      }

    input_stream >> keyword;
    DEBUG_ASSERT(keyword == "omega", "Expected keyword omega", keyword);

    std::getline(input_stream >> std::ws, line);

    std::istringstream omega_iss(line);
    values.clear();

    while (omega_iss >> value)
      {
        values.push_back(value);
      }

    if (values.size() == 1)
      {
        std::fill(omega.begin(), omega.end(), values[0]);
      }
    else
      {
        DEBUG_ASSERT(values.size() == N, "Expected 1 or N omega values", values.size());
        omega = values;
      }
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
