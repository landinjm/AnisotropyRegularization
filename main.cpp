#include <filesystem>
#include <fstream>
#include <iostream>

/**
 * Heaviside step function
 */
template <typename RealType>
inline constexpr RealType
theta(RealType x)
{
  return x >= RealType {0} ? RealType {1} : RealType {0};
}

template <typename RealType>
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
    alpha.resize(N);
    omega.resize(N);
  }

  /**
   * Validate the parameters
   */
  void
  validate() const
  {}
};

int
main(int argc, char *argv[])
{
  std::string   filename = (argc > 1) ? argv[1] : "parameter.prm";
  std::ifstream prm(filename);
  if (!prm)
    {
      std::cerr << "Could not open " << filename << '\n';
      return 1;
    }

  return 0;
}
