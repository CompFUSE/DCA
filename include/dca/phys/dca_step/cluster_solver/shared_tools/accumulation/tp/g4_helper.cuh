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

#include <array>
#include <cuda.h>
#include <vector>
#include <stdexcept>

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

class G4Helper {
public:
  __host__ void set(int nb, int nk, int nw_pos, const std::vector<int>& k_ex_indices,
                    const std::vector<int>& w_ex_indices, const int* add_k, int lda,
                    const int* sub_k, int lds);

  G4Helper(const G4Helper& other) = default;

  __host__ bool isInitialized() const {
    return device_members_ != nullptr;
  }

  // Returns the index of k + k_ex.
  __device__ inline int addKex(int k_idx, int k_ex_idx) const;
  // Returns the index of k_ex - k.
  __device__ inline int kexMinus(int k_idx, int k_ex_idx) const;
  // Returns the index of w + w_ex.
  __device__ inline int addWex(int w_idx, int w_ex_idx) const;
  // Returns the index of w_ex - w.
  __device__ inline int wexMinus(int w_idx, int w_ex_idx) const;

  // Maps the indices w1 w2 from the compact frequency domain of G4,
  // to the extended (positive for w1) domain used by G.
  // In/Out: w1, w2.
  // Returns: true if G(w1, w2) is stored as a complex conjugate.
  __device__ inline bool extendWIndices(int& w1, int& w2) const;

  // Returns the linear index of G4 as a function of
  // band, band, band, band, k1, k2, k_ex, w1, w2, w_ex.
  __device__ inline int g4Index(int b1, int b2, int b3, int b4, int k1, int k2, int k_ex, int w1,
                                int w2, int w_ex) const;
  // Single band version of the above method.
  __device__ inline int g4Index(int k1, int k2, int k_ex, int w1, int w2, int w_ex) const;

protected:
  // This object can be constructed only through its derived class.
  G4Helper() = default;

  int* add_matrix_ = nullptr;
  int* sub_matrix_ = nullptr;

  // device_members_ points to an array of private members stored in a single GPU allocation.
  // The order is:
  // lda: leading dimension of add_matrix_.
  // lds: leading dimension of sub_matrix_.
  // nw_pos: number of positive frequencies stored in G4.
  // ext_size: difference between the number of positive frequencies stored in G and G4.
  int* device_members_ = nullptr;

  int* w_ex_indices_ = nullptr;
  int* k_ex_indices_ = nullptr;

  // Stores the steps between each subdomain used by g4Index.
  int* sbdm_steps_ = nullptr;
};

class G4HelperManager : public G4Helper {
public:
  G4HelperManager() = default;
  __host__ ~G4HelperManager();
  __host__ __device__ G4HelperManager(const G4HelperManager& other) = delete;
};

__host__ void G4Helper::set(int nb, int nk, int nw_pos, const std::vector<int>& delta_k,
                            const std::vector<int>& delta_w, const int* add_k, int lda,
                            const int* sub_k, int lds) {
  if (isInitialized())
    throw(std::logic_error("already initialized."));

  cudaMalloc(&add_matrix_, sizeof(int) * lda * nk);
  cudaMalloc(&sub_matrix_, sizeof(int) * lds * nk);
  cudaMemcpy(add_matrix_, add_k, sizeof(int) * lda * nk, cudaMemcpyHostToDevice);
  cudaMemcpy(sub_matrix_, sub_k, sizeof(int) * lds * nk, cudaMemcpyHostToDevice);

  int ext_size = 0;
  for (const int idx : delta_w)
    ext_size = std::max(ext_size, std::abs(idx));

  const std::array<int, 4> device_members_host{lda, lds, nw_pos, ext_size};
  cudaMalloc(&device_members_, sizeof(int) * device_members_host.size());
  cudaMemcpy(device_members_, device_members_host.data(), sizeof(int) * device_members_host.size(),
             cudaMemcpyHostToDevice);

  cudaMalloc(&w_ex_indices_, sizeof(int) * delta_w.size());
  cudaMemcpy(w_ex_indices_, delta_w.data(), sizeof(int) * delta_w.size(), cudaMemcpyHostToDevice);
  cudaMalloc(&k_ex_indices_, sizeof(int) * delta_k.size());
  cudaMemcpy(k_ex_indices_, delta_k.data(), sizeof(int) * delta_k.size(), cudaMemcpyHostToDevice);

  const int nb4 = nb * nb * nb * nb;
  const int nk3 = nk * nk * delta_k.size();
  const int nw = 2 * nw_pos;
  const std::array<int, 10> steps_host{1,
                                       nb,
                                       nb * nb,
                                       nb * nb * nb,
                                       nb4,
                                       nb4 * nk,
                                       nb4 * nk * nk,
                                       nb4 * nk3,
                                       nb4 * nk3 * nw,
                                       nb4 * nk3 * nw * nw};
  cudaMalloc(&sbdm_steps_, sizeof(int) * steps_host.size());
  cudaMemcpy(sbdm_steps_, steps_host.data(), sizeof(int) * steps_host.size(), cudaMemcpyHostToDevice);
}

__device__ int G4Helper::addWex(const int w_idx, const int w_ex_idx) const {
  return w_idx + w_ex_indices_[w_ex_idx];
}

__device__ int G4Helper::wexMinus(const int w_idx, const int w_ex_idx) const {
  const int nw = 2 * device_members_[2];
  return w_ex_indices_[w_ex_idx] + nw - 1 - w_idx;
}

__device__ int G4Helper::addKex(const int k_idx, const int k_ex_idx) const {
  const int ld = device_members_[0];
  const int k_ex = k_ex_indices_[k_ex_idx];
  return add_matrix_[k_idx + ld * k_ex];
}

__device__ int G4Helper::kexMinus(const int k_idx, const int k_ex_idx) const {
  const int ld = device_members_[1];
  const int k_ex = k_ex_indices_[k_ex_idx];
  return sub_matrix_[k_idx + ld * k_ex];
}

__device__ bool G4Helper::extendWIndices(int& w1, int& w2) const {
  const int extension = device_members_[3];
  const int n_w_ext_pos = extension + device_members_[2];
  w1 += extension;
  w2 += extension;
  if (w1 >= n_w_ext_pos) {
    w1 -= n_w_ext_pos;
    return false;
  }
  else {
    w1 = n_w_ext_pos - 1 - w1;
    w2 = 2 * n_w_ext_pos - 1 - w2;
    return true;
  }
}

__device__ int G4Helper::g4Index(int b1, int b2, int b3, int b4, int k1, int k2, int k_ex, int w1,
                                 int w2, int w_ex) const {
  return sbdm_steps_[0] * b1 + sbdm_steps_[1] * b2 + sbdm_steps_[2] * b3 + sbdm_steps_[3] * b4 +
         sbdm_steps_[4] * k1 + sbdm_steps_[5] * k2 + sbdm_steps_[6] * k_ex + sbdm_steps_[7] * w1 +
         sbdm_steps_[8] * w2 + sbdm_steps_[9] * w_ex;
  ;
}

__device__ int G4Helper::g4Index(int k1, int k2, int k_ex, int w1, int w2, int w_ex) const {
  return sbdm_steps_[4] * k1 + sbdm_steps_[5] * k2 + sbdm_steps_[6] * k_ex + sbdm_steps_[7] * w1 +
         sbdm_steps_[8] * w2 + sbdm_steps_[9] * w_ex;
}

__host__ G4HelperManager::~G4HelperManager() {
  cudaFree(G4Helper::add_matrix_);
  cudaFree(G4Helper::sub_matrix_);
  cudaFree(G4Helper::device_members_);
  cudaFree(G4Helper::sbdm_steps_);
  cudaFree(G4Helper::w_ex_indices_);
  cudaFree(G4Helper::k_ex_indices_);
}

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH
