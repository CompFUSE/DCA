// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// CT-AUX walker tools class.
// It is templated on the device type (CPU|GPU).

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_HPP

#include <cassert>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

// BLAS
#include "dca/linalg/blas/blas1.hpp"
#include "dca/linalg/blas/blas2.hpp"
#include "dca/linalg/blas/blas3.hpp"

#include "dca/linalg/lapack/bennet_update.hpp"
#include "dca/linalg/lapack/inverse.hpp"
#include "dca/linalg/lapack/lapack.hpp"
#include "dca/linalg/lapack/solve.hpp"

#ifdef DCA_HAVE_GPU
// CUBLAS
#include "dca/linalg/blas/cublas1.hpp"
#include "dca/linalg/blas/cublas3.hpp"
#include "dca/linalg/blas/cublas_conversion_char_types.hpp"
#include "dca/linalg/blas/kernels_gpu.hpp"

#include "dca/linalg/lapack/laset_gpu.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/lapack/multiply_diagonal_gpu.hpp"

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/ct_aux_walker_tools_kernels.hpp"
#endif

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

//
// Empty class template
//
template <dca::linalg::DeviceType device_t, typename Scalar>
class CT_AUX_WALKER_TOOLS {};

//
// Specialization for CPU
// INTERNAL: Should we inline the static methods?
//
template <typename Scalar>
class CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar> {
  using Real = dca::util::RealAlias<Scalar>;
  const static int BLOCK_SIZE = 32;  // dca::linalg::Matrix<double, dca::linalg::CPU>::BLOCK_SIZE;

public:
  CT_AUX_WALKER_TOOLS(int k_ph);

  static void compute_Gamma(dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma,
                            dca::linalg::Matrix<Scalar, dca::linalg::CPU>& N,
                            dca::linalg::Matrix<Scalar, dca::linalg::CPU>& G_precomputed,
                            dca::linalg::Vector<int, dca::linalg::CPU>& random_vertex_vector,
                            dca::linalg::Vector<Scalar, dca::linalg::CPU>& exp_V,
                            dca::linalg::Vector<Scalar, dca::linalg::CPU>& exp_delta_V,
                            int thread_id, int stream_id);

  static void set_to_identity(dca::linalg::Matrix<Scalar, dca::linalg::CPU>& M, int index);

  // inline Scalar solve_Gamma(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU,
  // Scalar
  // exp_delta_V);
  auto solve_Gamma(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU,
                   Scalar exp_delta_V, Real& max, Real& min) -> Real;

  /** Solve gamma using blocked submatrix updates
   *  \param[in] int
   *  \param[inout] Gamma_LU
   *  \param[ioout] max
   *  \param[inout] min
   */
  auto solve_Gamma_blocked(const int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU,
                           Scalar exp_delta_V, Real& max, Real& min) -> Scalar;

  auto apply_bennett_on_Gamma(int k, int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU,
                              Scalar phani_gamma, Real& max, Real& min) -> Real;

  /** Debugging function
   *  requires: n > 0
   */
  static bool test_max_min(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU, Real max,
                    Real min);

private:
  void solve_Gamma_slow(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU);
  void solve_Gamma_slow(int n, Scalar* Gamma_LU, int lda);
  void solve_Gamma_fast(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU);
  void solve_Gamma_fast(int n, Scalar* A, int LD);
  void solve_Gamma_BLAS(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU);
  void solve_Gamma_BLAS(int n, Scalar* Gamma_LU, int lda);


  void solve_Gamma_blocked(int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU);

  /** Return x/y.
   *  nothing but that is done for real number but for complex division an implementation consistent with cuda
   *  and magma libraries is used.
   *  This is more overflow resistant than the std lib method.
   *  But we need this for reasonable numerical agreement for testing.
   */
  static Scalar consistentScalarDiv(Scalar x, Scalar y);

private:
  dca::linalg::Vector<Scalar, dca::linalg::CPU> r;
  dca::linalg::Vector<Scalar, dca::linalg::CPU> c;
  dca::linalg::Vector<Scalar, dca::linalg::CPU> d;
};

#ifdef DCA_HAVE_GPU
//
// Specialization for GPU
//
template <typename Scalar>
class CT_AUX_WALKER_TOOLS<dca::linalg::GPU, Scalar> {
public:
  static void compute_Gamma(dca::linalg::Matrix<Scalar, dca::linalg::GPU>& Gamma,
                            dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N,
                            dca::linalg::Matrix<Scalar, dca::linalg::GPU>& G,
                            dca::linalg::Vector<int, dca::linalg::GPU>& random_vertex_vector,
                            dca::linalg::Vector<Scalar, dca::linalg::GPU>& exp_V,
                            dca::linalg::Vector<Scalar, dca::linalg::GPU>& exp_delta_V,
                            int thread_id, int stream_id) {
    Gamma.resize(random_vertex_vector.size());

    assert(Gamma.nrRows() == Gamma.nrCols());

    walkerkernels::compute_Gamma(Gamma.ptr(), Gamma.nrRows(), Gamma.leadingDimension(), N.ptr(),
                                 N.nrRows(), N.nrCols(), N.leadingDimension(), G.ptr(), G.nrRows(),
                                 G.nrCols(), G.leadingDimension(), random_vertex_vector.ptr(),
                                 exp_V.ptr(), exp_delta_V.ptr(), thread_id, stream_id);
  }
};

#endif  // DCA_HAVE_GPU

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_HPP
