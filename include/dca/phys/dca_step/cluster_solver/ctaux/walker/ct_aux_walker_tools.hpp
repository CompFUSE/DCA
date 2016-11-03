// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// CT-AUX walker tools class.
// It is templated on the device type (CPU|GPU).

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_HPP

#include <cassert>

#include "dca/linalg/linalg.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

//
// Empty class template
//
template <dca::linalg::DeviceType device_t>
class CT_AUX_WALKER_TOOLS {};

//
// Specialization for CPU
// INTERNAL: Should we inline the static methods?
//
template <>
class CT_AUX_WALKER_TOOLS<dca::linalg::CPU> {
  const static int BLOCK_SIZE = 32;  // dca::linalg::Matrix<double, dca::linalg::CPU>::BLOCK_SIZE;

public:
  CT_AUX_WALKER_TOOLS(int k_ph);

  static void compute_Gamma(dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma,
                            dca::linalg::Matrix<double, dca::linalg::CPU>& N,
                            dca::linalg::Matrix<double, dca::linalg::CPU>& G_precomputed,
                            dca::linalg::Vector<int, dca::linalg::CPU>& random_vertex_vector,
                            dca::linalg::Vector<double, dca::linalg::CPU>& exp_V,
                            dca::linalg::Vector<double, dca::linalg::CPU>& exp_delta_V,
                            int thread_id, int stream_id);

  static void set_to_identity(dca::linalg::Matrix<double, dca::linalg::CPU>& M, int index);

  // inline double solve_Gamma(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU,
  // double
  // exp_delta_V);
  double solve_Gamma(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU,
                     double exp_delta_V, double& max, double& min);
  double solve_Gamma_blocked(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU,
                             double exp_delta_V, double& max, double& min);

  double apply_bennett_on_Gamma(int k, int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU,
                                double phani_gamma, double& max, double& min);

private:
  void solve_Gamma_slow(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU);
  void solve_Gamma_fast(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU);
  void solve_Gamma_BLAS(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU);

  void solve_Gamma_fast(int n, double* A, int LD);

  void solve_Gamma_blocked(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU);

  bool test_max_min(int n, dca::linalg::Matrix<double, dca::linalg::CPU>& Gamma_LU, double max,
                    double min);

private:
  dca::linalg::Vector<double, dca::linalg::CPU> r;
  dca::linalg::Vector<double, dca::linalg::CPU> c;
  dca::linalg::Vector<double, dca::linalg::CPU> d;
};

#ifdef DCA_HAVE_CUDA
//
// Specialization for GPU
//
namespace walkerkernels {
// dca::phys::solver::ctaux::walkerkernels::

void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r, int N_c, int N_ld,
                   double* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, double* exp_V,
                   double* exp_delta_V, int thread_id, int stream_id);

}  // walkerkernels
// dca::phys::solver::ctaux::

template <>
class CT_AUX_WALKER_TOOLS<dca::linalg::GPU> {
public:
  static void compute_Gamma(dca::linalg::Matrix<double, dca::linalg::GPU>& Gamma,
                            dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                            dca::linalg::Matrix<double, dca::linalg::GPU>& G,
                            dca::linalg::Vector<int, dca::linalg::GPU>& random_vertex_vector,
                            dca::linalg::Vector<double, dca::linalg::GPU>& exp_V,
                            dca::linalg::Vector<double, dca::linalg::GPU>& exp_delta_V,
                            int thread_id, int stream_id) {
    Gamma.resize(random_vertex_vector.size());

    assert(Gamma.nrRows() == Gamma.nrCols());

    walkerkernels::compute_Gamma(Gamma.ptr(), Gamma.nrRows(), Gamma.leadingDimension(), N.ptr(),
                                 N.nrRows(), N.nrCols(), N.leadingDimension(), G.ptr(), G.nrRows(),
                                 G.nrCols(), G.leadingDimension(), random_vertex_vector.ptr(),
                                 exp_V.ptr(), exp_delta_V.ptr(), thread_id, stream_id);
  }
};

#endif  // DCA_HAVE_CUDA

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_HPP
