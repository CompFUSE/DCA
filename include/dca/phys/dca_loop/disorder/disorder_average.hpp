// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira (anirudha.mirmira@gmail.com)
//
// Numeric kernels of the DISORDERED_G0 disorder-averaging path, factored out of the cluster
// solver and the DCA loop so they can be unit tested in isolation:
//
//   accumulateDisorderDyson - per-frequency dense (nu*N_R) Dyson G = G0_dis - G0_dis*M*G0_dis/beta,
//                             accumulated weighted into a two space-index function. Used by
//                             CtauxClusterSolver::accumulateGrrwFromMrrw.
//   translationalAverage    - project the disorder-averaged two-site G(<nu,R,nu,R,w>) onto the
//                             single-displacement g(<nu,nu,R,w>) with the cluster subtract table.
//                             Used by DcaLoop::averageGkw.

#ifndef DCA_PHYS_DCA_LOOP_DISORDER_DISORDER_AVERAGE_HPP
#define DCA_PHYS_DCA_LOOP_DISORDER_DISORDER_AVERAGE_HPP

#include <complex>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/linalg/matrixop.hpp"

namespace dca::phys::disorder {

// Per-frequency disorder Dyson on the dense two-site basis. At each w the (nu,R) indices form a
// (nu*N_R) square block (row index nu + nu_size*R, leading dimension nu*N_R; the contiguous layout
// of <nu,R,nu,R,w> functions). Computes G = G0_dis - G0_dis*M*G0_dis/beta and accumulates
// weight*G into G_acc (G_acc is NOT zeroed - the caller accumulates over configurations).
//
// G0_dis and m are inputs (taken by non-const ref so MatrixView can alias their storage); G_acc is
// the running accumulator. FG0/FM/FAcc are the (possibly distinct) two space-index function types.
template <class NuDmn, class RDmn, class WDmn, typename Real, class FG0, class FM, class FAcc>
void accumulateDisorderDyson(FG0& G0_dis, FM& m, double beta, double weight, FAcc& G_acc) {
  const std::size_t matrix_size = NuDmn::dmn_size() * RDmn::dmn_size();
  dca::linalg::Matrix<std::complex<Real>, dca::linalg::CPU> G0_times_M(matrix_size);
  dca::linalg::Matrix<std::complex<Real>, dca::linalg::CPU> G0_M_G0(matrix_size);
  using MatrixView = dca::linalg::MatrixView<std::complex<Real>, dca::linalg::CPU>;

  for (int w = 0; w < WDmn::dmn_size(); ++w) {
    const MatrixView G0_matrix(&G0_dis(0, 0, 0, 0, w), matrix_size);
    const MatrixView M_matrix(&m(0, 0, 0, 0, w), matrix_size);

    dca::linalg::matrixop::gemm(G0_matrix, M_matrix, G0_times_M);  // G0 * M
    dca::linalg::matrixop::gemm(G0_times_M, G0_matrix, G0_M_G0);   // (G0 * M) * G0

    MatrixView acc(&G_acc(0, 0, 0, 0, w), matrix_size);
    for (int j = 0; j < static_cast<int>(matrix_size); ++j)
      for (int i = 0; i < static_cast<int>(matrix_size); ++i)
        acc(i, j) += static_cast<Real>(weight) * (G0_matrix(i, j) - G0_M_G0(i, j) / beta);
  }
}

// Translational average of the disorder-averaged two-site G: project <nu1,R1,nu2,R2,w> onto the
// single-displacement g <nu1,nu2,R,w>. Disorder averaging restores translational symmetry, so g
// depends only on the displacement R. Summing over R1 = R0 at fixed displacement R uses
// R2 = subtract(R,R0) (subtract(i,j) = index of r_j - r_i), the exact inverse of the unfold
// g(R1-R2) -> G(R1,R2) that uses subtract(R2,R1). g_r is overwritten (zeroed then filled).
template <class NuDmn, class RDmn, class WDmn, class TwoRFunction, class OneRFunction>
void translationalAverage(const TwoRFunction& g_rr, OneRFunction& g_r) {
  const int N = RDmn::dmn_size();
  g_r = 0;
  for (int w = 0; w < WDmn::dmn_size(); ++w)
    for (int R = 0; R < N; ++R)
      for (int n2 = 0; n2 < NuDmn::dmn_size(); ++n2)
        for (int n1 = 0; n1 < NuDmn::dmn_size(); ++n1)
          for (int R0 = 0; R0 < N; ++R0)
            g_r(n1, n2, R, w) += g_rr(n1, R0, n2, RDmn::parameter_type::subtract(R, R0), w);
  g_r /= static_cast<double>(N);
}

}  // namespace dca::phys::disorder

#endif  // DCA_PHYS_DCA_LOOP_DISORDER_DISORDER_AVERAGE_HPP
