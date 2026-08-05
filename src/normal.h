#pragma once

#include <cmath>
#include <libassert/assert.hpp>

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
  {
    using std::cos;
    using std::sin;

    x = cos(theta);
    y = sin(theta);
  }

  Normal(RealType x_, RealType y_)
  {
    using std::sqrt;

    RealType d = sqrt(x_ * x_ + y_ * y_);
    DEBUG_ASSERT(d > RealType {0}, "Vector input length must be nonzero");

    x = x_ / d;
    y = y_ / d;
  }

  template <typename number>
  [[nodiscard]] constexpr RealType
  dot(const Normal<2, number> &other) const
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
  {
    using std::cos;
    using std::sin;

    x = cos(psi) * cos(theta);
    y = cos(psi) * sin(theta);
    z = sin(psi);
  }

  Normal(RealType x_, RealType y_, RealType z_)
  {
    using std::sqrt;

    RealType d = sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    DEBUG_ASSERT(d > RealType {0}, "Vector input length must be nonzero");

    x = x_ / d;
    y = y_ / d;
    z = z_ / d;
  }

  template <typename number>
  [[nodiscard]] constexpr RealType
  dot(const Normal<3, number> &other) const
  {
    return x * other.x + y * other.y + z * other.z;
  }
};
