#pragma once

#include "normal.h"

#include <fstream>
#include <libassert/assert.hpp>
#include <numbers>
#include <string>
#include <vector>

template <int dim, typename RealType, typename Functor>
void
write_gamma_surface_vtk(const std::string &filename,
                        const Functor     &func,
                        int                n_theta = 181,
                        int                n_psi   = 91)
{
  std::ofstream vtk(filename);
  DEBUG_ASSERT(vtk, "Could not open output file.");

  vtk << "# vtk DataFile Version 3.0\n";
  vtk << "Surface energy\n";
  vtk << "ASCII\n";
  vtk << "DATASET STRUCTURED_GRID\n";

  std::vector<RealType> gamma_values;

  if constexpr (dim == 2)
    {
      vtk << "DIMENSIONS " << n_theta << " 1 1\n";
      vtk << "POINTS " << n_theta << " double\n";

      gamma_values.reserve(n_theta);
      for (int i = 0; i < n_theta; ++i)
        {
          RealType              theta = 2.0 * std::numbers::pi * i / (n_theta - 1);
          Normal<dim, RealType> n(theta);
          RealType              gamma = func(n);
          gamma_values.push_back(gamma);
          vtk << gamma * n.x << ' ' << gamma * n.y << ' ' << " 0\n";
        }
      vtk << "\nPOINT_DATA " << n_theta << '\n';
    }
  else if constexpr (dim == 3)
    {
      vtk << "DIMENSIONS " << n_theta << " " << n_psi << " 1\n";
      vtk << "POINTS " << n_theta * n_psi << " double\n";

      gamma_values.reserve(n_theta * n_psi);
      for (int j = 0; j < n_psi; ++j)
        {
          RealType psi = -std::numbers::pi / 2.0 + std::numbers::pi * j / (n_psi - 1);
          for (int i = 0; i < n_theta; ++i)
            {
              RealType              theta = 2.0 * std::numbers::pi * i / (n_theta - 1);
              Normal<dim, RealType> n(theta, psi);
              RealType              gamma = func(n);
              gamma_values.push_back(gamma);
              vtk << gamma * n.x << ' ' << gamma * n.y << ' ' << gamma * n.z << '\n';
            }
        }
      vtk << "\nPOINT_DATA " << n_theta * n_psi << '\n';
    }
  else
    {
      static_assert(dim == 2 || dim == 3,
                    "Anisotropies are only supported for 2D and 3D");
    }

  vtk << "SCALARS gamma double\n";
  vtk << "LOOKUP_TABLE default\n";
  for (RealType g : gamma_values)
    {
      vtk << g << '\n';
    }
}

template <int dim, typename RealType, typename Functor>
void
write_wulff_surface_vtk(const std::string &filename,
                        const Functor     &func,
                        int                n_theta = 181,
                        int                n_psi   = 91)
{
  std::ofstream vtk(filename);
  DEBUG_ASSERT(vtk, "Could not open output file.");

  vtk << "# vtk DataFile Version 3.0\n";
  vtk << "Surface energy\n";
  vtk << "ASCII\n";
  vtk << "DATASET STRUCTURED_GRID\n";

  std::vector<RealType> gamma_values;

  if constexpr (dim == 2)
    {
      vtk << "DIMENSIONS " << n_theta << " 1 1\n";
      vtk << "POINTS " << n_theta << " double\n";

      gamma_values.reserve(n_theta);
      for (int i = 0; i < n_theta; ++i)
        {
          RealType              theta = 2.0 * std::numbers::pi * i / (n_theta - 1);
          Normal<dim, RealType> n(theta);
          Normal<dim, RealType> u     = func.wulff_point(n);
          RealType              gamma = func(n);
          gamma_values.push_back(gamma);
          vtk << u.x << ' ' << u.y << ' ' << " 0\n";
        }
      vtk << "\nPOINT_DATA " << n_theta << '\n';
    }
  else if constexpr (dim == 3)
    {
      vtk << "DIMENSIONS " << n_theta << " " << n_psi << " 1\n";
      vtk << "POINTS " << n_theta * n_psi << " double\n";

      gamma_values.reserve(n_theta * n_psi);
      for (int j = 0; j < n_psi; ++j)
        {
          RealType psi = -std::numbers::pi / 2.0 + std::numbers::pi * j / (n_psi - 1);
          for (int i = 0; i < n_theta; ++i)
            {
              RealType              theta = 2.0 * std::numbers::pi * i / (n_theta - 1);
              Normal<dim, RealType> n(theta, psi);
              Normal<dim, RealType> u     = func.wulff_point(n);
              RealType              gamma = func(n);
              gamma_values.push_back(gamma);
              vtk << u.x << ' ' << u.y << ' ' << u.z << '\n';
            }
        }
      vtk << "\nPOINT_DATA " << n_theta * n_psi << '\n';
    }
  else
    {
      static_assert(dim == 2 || dim == 3,
                    "Anisotropies are only supported for 2D and 3D");
    }

  vtk << "SCALARS gamma double\n";
  vtk << "LOOKUP_TABLE default\n";
  for (RealType g : gamma_values)
    {
      vtk << g << '\n';
    }
}
