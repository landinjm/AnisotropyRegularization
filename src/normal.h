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

  [[nodiscard]] constexpr RealType
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

  [[nodiscard]] constexpr RealType
  dot(const Normal &other) const
  {
    return x * other.x + y * other.y + z * other.z;
  }
};
