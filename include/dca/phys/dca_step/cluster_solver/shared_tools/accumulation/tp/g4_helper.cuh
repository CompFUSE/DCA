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
  static void set(int nb, int nk, int nw_pos, const std::vector<int>& k_ex_indices,
                  const std::vector<int>& w_ex_indices, const int* add_k, int lda, const int* sub_k,
                  int lds, int k0);

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

  // Returns the linear index of G4 as a function of
  // band, band, band, band, k1, k2, k_ex, w1, w2, w_ex.
  __device__ inline unsigned int g4Index(int b1, int b2, int b3, int b4, int k1, int k2, int k_ex,
                                         int w1, int w2, int w_ex) const;
  // Single band version of the above method.
  __device__ inline unsigned int g4Index(int k1, int k2, int k_ex, int w1, int w2, int w_ex) const;

  // Returns range (start and end index) of G4 in which local rank should compute, when GPUDirect is enabled
  __device__ inline void getComputeRange(int my_rank, int mpi_size, int total_G4_size, int& start, int& end) const;

protected:
  int lda_;
  int lds_;
  int nw_pos_;
  int ext_size_;
  unsigned sbdm_steps_[10];
  const int* w_ex_indices_;
  const int* k_ex_indices_;
};

// Global instance to be used in the tp accumulation kernel.
extern __device__ __constant__ G4Helper g4_helper;

inline __device__ int G4Helper::addWex(const int w_idx, const int w_ex_idx) const {
  return w_idx + w_ex_indices_[w_ex_idx];
}

inline __device__ int G4Helper::wexMinus(const int w_idx, const int w_ex_idx) const {
  const int nw = 2 * nw_pos_;
  return w_ex_indices_[w_ex_idx] + nw - 1 - w_idx;
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

inline __device__ unsigned int G4Helper::g4Index(int b1, int b2, int b3, int b4, int k1, int k2,
                                                 int k_ex, int w1, int w2, int w_ex) const {
  return sbdm_steps_[0] * b1 + sbdm_steps_[1] * b2 + sbdm_steps_[2] * b3 + sbdm_steps_[3] * b4 +
         sbdm_steps_[4] * k1 + sbdm_steps_[5] * k2 + sbdm_steps_[6] * k_ex + sbdm_steps_[7] * w1 +
         sbdm_steps_[8] * w2 + sbdm_steps_[9] * w_ex;
  ;
}

inline __device__ unsigned int G4Helper::g4Index(int k1, int k2, int k_ex, int w1, int w2,
                                                 int w_ex) const {
  return sbdm_steps_[4] * k1 + sbdm_steps_[5] * k2 + sbdm_steps_[6] * k_ex + sbdm_steps_[7] * w1 +
         sbdm_steps_[8] * w2 + sbdm_steps_[9] * w_ex;
}

#ifdef DCA_WITH_NVLINK
inline __device__
void G4Helper::getComputeRange(int my_rank, int mpi_size, int total_G4_size, int& start, int& end) const {

    unsigned int offset = 0;
    // check if originally flattened one-dimensional G4 array can be equally (up to 0) distributed across ranks
    // if balanced, each rank has same amount of elements to compute
    // if not, ranks with (rank_id < nb_more_work_ranks) has to compute 1 more element than other ranks
    bool balanced = (total_G4_size % mpi_size == 0);
    int local_work = total_G4_size / mpi_size;

    if(balanced) {
        offset = my_rank  * local_work;
        end  = offset + local_work;
    }
    else {
        int nb_more_work_ranks = total_G4_size % mpi_size;

        if (my_rank < nb_more_work_ranks) {
            offset = my_rank * (local_work + 1);
            end = offset + (local_work + 1);
        }
        else {
            offset = nb_more_work_ranks * (local_work + 1) + (my_rank - nb_more_work_ranks) * local_work;
            end = offset + local_work;
        }
    }
    start = offset;
}
#endif

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH
