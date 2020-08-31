// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Helper class for adding and subtracting momentum and frequency on the device.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH

#include <cassert>
#include <vector>

#include <cuda.h>

#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

class G4Helper {
public:
  static void set(unsigned int nb, unsigned int nk, unsigned int nw_pos,
                  const std::vector<int>& delta_k, const std::vector<int>& delta_w,
                  const int* add_k, unsigned int lda, const int* sub_k, unsigned int lds);

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
  // In/Out: w1, w2.
  // Returns: true if G(w1, w2) is stored as a complex conjugate.
  __device__ inline bool extendGIndices(int& k1, int& k2, int& w1, int& w2) const;

  // Unroll the linear index of G4 as a function of band, band, band, band,
  // k1, k2, k_ex, w1, w2, w_ex.
  __device__ inline void unrollIndex(std::size_t index, unsigned& b1, unsigned& b2, unsigned& b3,
                                     unsigned& b4, unsigned& k1, unsigned& w1, unsigned& k2,
                                     unsigned& w2, unsigned& k_ex, unsigned& w_ex) const;

protected:
  std::size_t sbdm_steps_[10];

  const int* w_ex_indices_;
  const int* k_ex_indices_;
  unsigned ext_size_;

  unsigned nw_pos_;
  unsigned nb_;
  unsigned nc_;
  unsigned nw_;
  unsigned n_k_ex_;
  unsigned n_w_ex_;
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

inline __device__ bool G4Helper::extendGIndices(int& k1, int& k2, int& w1, int& w2) const {
  const int n_w_ext_pos = ext_size_ + nw_pos_;
  w1 += ext_size_;
  w2 += ext_size_;
  if (w1 >= n_w_ext_pos) {
    w1 -= n_w_ext_pos;
    return false;
  }
  else {
    w1 = n_w_ext_pos - 1 - w1;
    w2 = 2 * n_w_ext_pos - 1 - w2;
    k1 = kMinus(k1);
    k2 = kMinus(k2);
    return true;
  }
}

__device__ inline void G4Helper::unrollIndex(std::size_t index, unsigned& b1, unsigned& b2,
                                             unsigned& b3, unsigned& b4, unsigned& k1, unsigned& w1,
                                             unsigned& k2, unsigned& w2, unsigned& k_ex,
                                             unsigned& w_ex) const {
  auto unroll = [&](const unsigned dimension) {
    unsigned result = index / sbdm_steps_[dimension];
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
