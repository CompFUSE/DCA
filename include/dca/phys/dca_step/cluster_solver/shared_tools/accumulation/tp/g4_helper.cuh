// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Weile Wei (wwei9@lsu.edu)
//         Peter Doak (doakpw@ornl.gov)
//
// Helper class for adding and subtracting momentum and frequency on the device.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH

#include <cassert>
#include <vector>

#if defined(DCA_HAVE_GPU)
#include "dca/platform/dca_gpu.h"
#endif
#include "dca/distribution/dist_types.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

class G4Helper {
public:
  /** pass a bunch of information into the G4Helper so it can help./
   *  \param[in] extension_offset    (WTpExtDmn::dmn_size() - WTpDmn::dmn_size()) / 2
   */
  static void set(int nb, int nk, int nw,
                  const std::vector<int>& delta_k, const std::vector<int>& delta_w,
                  const int extension_offset, const int* add_k, int lda, const int* sub_k, int lds);

  __device__ auto get_bands() const {
    return nb_;
  }
  __device__ auto get_cluster_size() const {
    return nc_;
  }

  // Returns the index of k + k_ex.
  __device__ inline int addKex(int k_idx, int k_ex_idx) const;
  // Returns the index of k_ex - k.
  __device__ inline int kexMinus(int k_idx, int k_ex_idx) const;
  // Returns the index of -k.
  __device__ inline int kMinus(int k_idx) const;
  // Returns the index of w + w_ex.
  __device__ inline int addWex(int w_idx, int w_ex_idx) const;
  // Returns the index of w_ex - w.
  __device__ inline int wexMinus(int w_idx, int w_ex_idx) const;

  // Maps the indices w1 w2 from the compact frequency domain of G4,
  // to the extended (positive for w1) domain used by G.
  // k1, k2 mapped to minusK(k1), minus(k2)
  // In/Out: k1, k2, w1, w2.
  // Returns: true if G(w1, w2) is stored as a complex conjugate.
  __device__ bool extendGIndices(int& k1, int& k2, int& w1, int& w2) const;

  // Maps the indices w1 w2 from the compact frequency domain of G4,
  // to the extended (positive for w1) domain used by G.
  // In/Out: k1, k2, w1, w2.
  // Returns: true if G(w1, w2) is stored as a complex conjugate.
  __device__ bool extendGIndicesMultiBand(int& k1, int& k2, int& w1, int& w2) const;

  // Unroll the linear index of G4 as a function of band, band, band, band,
  // k1, k2, k_ex, w1, w2, w_ex.
  __device__ inline void unrollIndex(std::size_t index, int& b1, int& b2, int& b3,
                                     int& b4, int& k1, int& w1, int& k2,
                                     int& w2, int& k_ex, int& w_ex) const;

protected:
  std::size_t sbdm_steps_[10];

  const int* w_ex_indices_;
  const int* k_ex_indices_;
  int ext_size_;

  int nw_;
  int nb_;
  int nc_;
  int n_k_ex_;
  int n_w_ex_;
#ifndef NDEBUG
  int* bad_indicies_;
#endif
};

// Global instance to be used in the tp accumulation kernel.
extern __device__ __constant__ G4Helper g4_helper;

inline __device__ int G4Helper::addWex(const int w_idx, const int w_ex_idx) const {
  return w_idx + w_ex_indices_[w_ex_idx];
}

inline __device__ int G4Helper::wexMinus(const int w_idx, const int w_ex_idx) const {
  return w_ex_indices_[w_ex_idx] + nw_ - 1 - w_idx;
}

inline __device__ int G4Helper::addKex(const int k_idx, const int k_ex_idx) const {
  const int k_ex = k_ex_indices_[k_ex_idx];
  return solver::details::cluster_momentum_helper.add(k_idx, k_ex);
}

inline __device__ int G4Helper::kexMinus(const int k_idx, const int k_ex_idx) const {
  const int k_ex = k_ex_indices_[k_ex_idx];
  return solver::details::cluster_momentum_helper.subtract(k_idx, k_ex);
}

inline __device__ int G4Helper::kMinus(const int k_idx) const {
  return solver::details::cluster_momentum_helper.minus(k_idx);
}


inline __device__ void G4Helper::unrollIndex(std::size_t index, int& b1, int& b2,
                                             int& b3, int& b4, int& k1, int& w1,
                                             int& k2, int& w2, int& k_ex,
                                             int& w_ex) const {
  auto unroll = [&](const int dimension) {
    int result = index / sbdm_steps_[dimension];
    index -= result * sbdm_steps_[dimension];
    return result;
  };

  w_ex = unroll(9);
  k_ex = unroll(8);
  w2 = unroll(7);
  k2 = unroll(6);
  w1 = unroll(5);
  k1 = unroll(4);
  b4 = unroll(3);
  b3 = unroll(2);
  b2 = unroll(1);
  b1 = unroll(0);
}

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH
