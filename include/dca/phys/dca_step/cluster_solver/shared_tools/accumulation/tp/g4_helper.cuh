// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
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
  __host__ void set(int nb, int nk, int nw_pos, const std::vector<int>& q_indices,
                    const std::vector<int>& dw_indices, const int* add_k, int lda, const int* sub_k,
                    int lds);

  G4Helper(const G4Helper& other) = default;

  __host__ bool isInitialized() const {
    return parameters_ != nullptr;
  }

  __device__ inline int addQ(int k_idx, int q_idx) const;
  __device__ inline int qMinus(int k_idx, int q_idx) const;
  __device__ inline int addW(int w_idx, int dw_idx) const;
  __device__ inline int wMinus(int w_idx, int dw_idx) const;

  __device__ inline bool extendWIndices(int& w1, int& w2) const;

  __device__ inline int g4Index(int b1, int b2, int b3, int b4, int k1, int k2, int k3, int w1,
                                int w2, int w3) const;
  __device__ inline int g4Index(int k1, int k2, int k3, int w1, int w2, int w3) const;

protected:
  // This object can be constructed only through its derived class.
  G4Helper() = default;

  int* add_matrix_ = nullptr;
  int* sub_matrix_ = nullptr;

  // Stores in this order { lda, lds, nw_pos, ext_size}
  int* parameters_ = nullptr;

  int* dw_indices_ = nullptr;
  int* q_indices_ = nullptr;
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

  const std::array<int, 4> parameters_host{lda, lds, nw_pos, ext_size};
  cudaMalloc(&parameters_, sizeof(int) * parameters_host.size());
  cudaMemcpy(parameters_, parameters_host.data(), sizeof(int) * parameters_host.size(),
             cudaMemcpyHostToDevice);

  cudaMalloc(&dw_indices_, sizeof(int) * delta_w.size());
  cudaMemcpy(dw_indices_, delta_w.data(), sizeof(int) * delta_w.size(), cudaMemcpyHostToDevice);
  cudaMalloc(&q_indices_, sizeof(int) * delta_k.size());
  cudaMemcpy(q_indices_, delta_k.data(), sizeof(int) * delta_k.size(), cudaMemcpyHostToDevice);

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

__device__ int G4Helper::addW(const int w_idx, const int dw_idx) const {
  return w_idx + dw_indices_[dw_idx];
}

__device__ int G4Helper::wMinus(const int w_idx, const int dw_idx) const {
  const int nw = 2 * parameters_[2];
  return dw_indices_[dw_idx] + nw - 1 - w_idx;
}

__device__ int G4Helper::addQ(const int k_idx, const int q_idx) const {
  const int ld = parameters_[0];
  const int q = q_indices_[q_idx];
  return add_matrix_[k_idx + ld * q];
}

__device__ int G4Helper::qMinus(const int k_idx, const int q_idx) const {
  const int ld = parameters_[1];
  const int q = q_indices_[q_idx];
  return sub_matrix_[q + ld * k_idx];
}

__device__ bool G4Helper::extendWIndices(int& w1, int& w2) const {
  const int extension = parameters_[3];
  const int n_w_ext_pos = extension + parameters_[2];
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

__device__ int G4Helper::g4Index(int b1, int b2, int b3, int b4, int k1, int k2, int kex, int w1,
                                 int w2, int wex) const {
  return sbdm_steps_[0] * b1 + sbdm_steps_[1] * b2 + sbdm_steps_[2] * b3 + sbdm_steps_[3] * b4 +
         sbdm_steps_[4] * k1 + sbdm_steps_[5] * k2 + sbdm_steps_[6] * kex + sbdm_steps_[7] * w1 +
         sbdm_steps_[8] * w2 + sbdm_steps_[9] * wex;
  ;
}

__device__ int G4Helper::g4Index(int k1, int k2, int kex, int w1, int w2, int wex) const {
  return sbdm_steps_[4] * k1 + sbdm_steps_[5] * k2 + sbdm_steps_[6] * kex + sbdm_steps_[7] * w1 +
         sbdm_steps_[8] * w2 + sbdm_steps_[9] * wex;
}

__host__ G4HelperManager::~G4HelperManager() {
  cudaFree(G4Helper::add_matrix_);
  cudaFree(G4Helper::sub_matrix_);
  cudaFree(G4Helper::parameters_);
  cudaFree(G4Helper::sbdm_steps_);
  cudaFree(G4Helper::dw_indices_);
  cudaFree(G4Helper::q_indices_);
}

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_HELPER_CUH
