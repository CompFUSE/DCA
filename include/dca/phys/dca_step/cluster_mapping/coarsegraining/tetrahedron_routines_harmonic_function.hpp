// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Tetrahedron routines: harmonic function.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_ROUTINES_HARMONIC_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_ROUTINES_HARMONIC_FUNCTION_HPP

#include <complex>
#include <vector>

#include "dca/math/geometry/tetrahedron_mesh/tetrahedron.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

class tetrahedron_routines_harmonic_function {
public:
  // 1D
  static std::complex<double> execute(const std::vector<double>& r_vec,
                                      const math::geometry::tetrahedron<1>& tetrahedron);
  // 2D
  static std::complex<double> execute(const std::vector<double>& r_vec,
                                      const math::geometry::tetrahedron<2>& tetrahedron);

  // 3D
  static std::complex<double> execute(const std::vector<double>& r_vec,
                                      const math::geometry::tetrahedron<3>& tetrahedron);

private:
  // 2D cases:
  static void permute(math::geometry::tetrahedron<2>& tetrahedron_new,
                      const math::geometry::tetrahedron<2>& tetrahedron_old);

  static std::complex<double> case_2D(double dotRD1, double dotRD2, double dotRD2minD1);
  static std::complex<double> case_d1_2D(double dotRD1, double dotRD2, double dotRD2minD1);
  static std::complex<double> case_d2_2D(double dotRD1, double dotRD2, double dotRD2minD1);

  // 3D
  static void permute(math::geometry::tetrahedron<3>& tetrahedron_new,
                      const math::geometry::tetrahedron<3>& tetrahedron_old);

  static std::complex<double> case_3D(double dotRD1, double dotRD2, double dotRD3,
                                      double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

  static std::complex<double> case_d1_3D(double dotRD1, double dotRD2, double dotRD3,
                                         double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

  static std::complex<double> case_d2_3D(double dotRD1, double dotRD2, double dotRD3,
                                         double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

  static std::complex<double> case_d3_3D(double dotRD1, double dotRD2, double dotRD3,
                                         double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

  static std::complex<double> case_d1_d2_3D(double dotRD1, double dotRD2, double dotRD3,
                                            double dotRD2minD1, double dotRD3minD2,
                                            double dotRD1minD3);

  static std::complex<double> case_d2_d3_3D(double dotRD1, double dotRD2, double dotRD3,
                                            double dotRD2minD1, double dotRD3minD2,
                                            double dotRD1minD3);

  static std::complex<double> case_d3_d1_3D(double dotRD1, double dotRD2, double dotRD3,
                                            double dotRD2minD1, double dotRD3minD2,
                                            double dotRD1minD3);
};

}  // namespace clustermapping
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_ROUTINES_HARMONIC_FUNCTION_HPP
