#include <array>
#include <cmath>

template <uint dim, typename Real>
class Grad;

template <uint dim, typename Real>
Real
anisotropy(Grad<dim, Real> &grad);

template <uint dim, typename Real>
Real
grad_to_theta(Grad<dim, Real> &grad);

template <typename Real>
Real
grad_to_psi(Grad<3, Real> &grad);

/**
 * @brief User defined anisotropy function.
 *
 * NOTE: This function must be defined for a specific dimension and real type (e.g., float
 * or double).
 */
template <>
double
anisotropy<2, double>(Grad<2, double> &grad)
{
  auto theta = grad_to_theta(grad);

  return 1.0 + 0.1 * std::cos(4.0 * theta);
}

template <uint dim, typename Real>
class Grad
{
public:
  static_assert(dim == 2 || dim == 3, "Only 2D and 3D are supported");
  Grad() = default;

  template <typename... Args>
  explicit Grad(Args... args)
    : data {static_cast<Real>(args)...}
  {
    static_assert(sizeof...(Args) == dim, "Wrong number of arguments");
  }

  Real &
  operator[](std::size_t i)
  {
    return data[i];
  }

  const Real &
  operator[](std::size_t i) const
  {
    return data[i];
  }

private:
  std::array<Real, dim> data;
};

template <uint dim, typename Real>
Real
grad_to_theta(Grad<dim, Real> &grad)
{
  static_assert(dim == 2 || dim == 3, "Only 2D and 3D are supported");
  return std::atan2(grad[1], grad[0]);
}

template <typename Real>
Real
grad_to_psi(Grad<3, Real> &grad)
{
  return std::atan2(std::sqrt(grad[0] * grad[0] + grad[1] * grad[1]), grad[2]);
}

template <uint dim, typename Real>
find_missing_angles()
{
  static_assert(dim == 2 || dim == 3, "Only 2D and 3D are supported");

  if constexpr (dim == 2)
    {
    }
  else if constexpr (dim == 3)
    {
    }
}

int
main()
{
  find_missing_angles<2, double>();

  return 0;
}
