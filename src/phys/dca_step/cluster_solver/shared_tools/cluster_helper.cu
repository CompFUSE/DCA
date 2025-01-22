// Copyright (C) 2025 ETH Zurich
// Copyright (C) 2025 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This file implements ClusterHelper::set.

#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

#include <mutex>
#include <array>
#include <cassert>
#include "dca/platform/dca_gpu.h"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

__device__ __constant__ ClusterHelper cluster_real_helper;
__device__ __constant__ ClusterHelper cluster_momentum_helper;
__CONSTANT__ int* cluster_add_matrix;
__CONSTANT__ int* cluster_sub_matrix;
__CONSTANT__ int* cluster_momentum_add_matrix;
__CONSTANT__ int* cluster_momentum_sub_matrix;

__global__ void checkClusterHelper(int nc, int lds) {
  const int i_t = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i_t == 0 && j == 0) {
    printf("\nClusterHelpMatrix check:\n");
    for (int ch_i = 0; ch_i < nc; ++ch_i) {
      for (int ch_j = 0; ch_j < lds; ++ch_j) {
        printf("%d ", cluster_real_helper.sub_matrix_[ch_i * lds + ch_j]);
      }
      printf("\n");
    }
  }
}

__global__ void checkClusterMomentumHelper(int nc, int lds) {
  const int i_t = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i_t == 0 && j == 0) {
    printf("ClusterMomentumHelpMatrix check:\n");
    for (int ch_i = 0; ch_i < nc; ++ch_i) {
      for (int ch_j = 0; ch_j < lds; ++ch_j) {
        printf("%d ", cluster_momentum_helper.sub_matrix_[ch_i * lds + ch_j]);
      }
      printf("\n");
    }
  }
}

void ClusterHelper::setMomentum(int nc, const int* add, int lda, const int* sub, int lds) {
  static std::once_flag flag;

  std::call_once(flag, [=]() {
    ClusterHelper host_helper;
    host_helper.nc_ = nc;

    auto compact_transfer = [=](const int* matrix, int ldm, int** dest) {
      linalg::util::HostVector<int> compact(nc * nc);
      for (int j = 0; j < nc; ++j)
        for (int i = 0; i < nc; ++i)
          compact[i + nc * j] = matrix[i + ldm * j];

      cudaMalloc(dest, sizeof(int) * lds * nc);
      checkRC(cudaMemcpy(*dest, compact.data(), sizeof(int) * nc * nc, cudaMemcpyHostToDevice));
    };

    compact_transfer(add, lda, const_cast<int**>(&host_helper.add_matrix_));
    compact_transfer(sub, lds, const_cast<int**>(&host_helper.sub_matrix_));
    // The logic here for CUDA is transfers from unpinned memory only block for the copy from host
    // memory to a DMA buffer.  This is on the default queue here and so a stream synch will not
    // work to insure the ClusterHelper copy is actually complete.
    checkRC(cudaDeviceSynchronize());

    size_t cluster_helper_size;

    checkRC(cudaMemcpyToSymbol(HIP_SYMBOL(cluster_momentum_add_matrix), &host_helper.add_matrix_,
                               sizeof(int*)));
    checkRC(cudaMemcpyToSymbol(HIP_SYMBOL(cluster_momentum_sub_matrix), &host_helper.sub_matrix_,
                               sizeof(int*)));
    checkRC(cudaGetSymbolSize(&cluster_helper_size, HIP_SYMBOL(cluster_momentum_helper)));
    assert(cluster_helper_size == sizeof(ClusterHelper));
    checkRC(cudaMemcpyToSymbol(HIP_SYMBOL(cluster_momentum_helper), &host_helper,
                               sizeof(ClusterHelper)));
    checkRC(cudaDeviceSynchronize());
#ifdef DEBUG_CLUSTER_HELPER
    checkClusterMomentumHelper<<<1, 1, 0, 0>>>(nc, lds);
#endif
  });
}

void ClusterHelper::set(int nc, const int* add, int lda, const int* sub, int lds) {
  static std::once_flag flag;

  std::call_once(flag, [=]() {
    ClusterHelper host_helper;
    host_helper.nc_ = nc;

    auto compact_transfer = [=](const int* matrix, int ldm, int** dest) {
      linalg::util::HostVector<int> compact(nc * nc);
      for (int j = 0; j < nc; ++j)
        for (int i = 0; i < nc; ++i)
          compact[i + nc * j] = matrix[i + ldm * j];

      cudaMalloc(dest, sizeof(int) * lds * nc);
      checkRC(cudaMemcpy(*dest, compact.data(), sizeof(int) * nc * nc, cudaMemcpyHostToDevice));
    };

    compact_transfer(add, lda, const_cast<int**>(&host_helper.add_matrix_));
    compact_transfer(sub, lds, const_cast<int**>(&host_helper.sub_matrix_));
    // The logic here for CUDA is transfers from unpinned memory only block for the copy from host
    // memory to a DMA buffer.  This is on the default queue here and so a stream synch will not
    // work to insure the ClusterHelper copy is actually complete.
    checkRC(cudaDeviceSynchronize());

    size_t cluster_helper_size;

    // In debug on sdgx-2 for CTINT I see know evidence this actually works, it appears not to.
    checkRC(cudaGetSymbolSize(&cluster_helper_size, HIP_SYMBOL(cluster_real_helper)));
    assert(cluster_helper_size == sizeof(ClusterHelper));

    checkRC(
        cudaMemcpyToSymbol(HIP_SYMBOL(cluster_add_matrix), &host_helper.add_matrix_, sizeof(int*)));
    checkRC(
        cudaMemcpyToSymbol(HIP_SYMBOL(cluster_sub_matrix), &host_helper.sub_matrix_, sizeof(int*)));
    checkRC(cudaMemcpyToSymbol(HIP_SYMBOL(cluster_real_helper), &host_helper, cluster_helper_size));
    checkRC(cudaDeviceSynchronize());

#ifdef DEBUG_CLUSTER_HELPER
    checkClusterHelper<<<1, 1, 0, 0>>>(nc, lds);
#endif
    }
    else {
      // In debug on sdgx-2 for CTINT I see know evidence this actually works, it appears not to.
      size_t cluster_helper_size;
      checkRC(cudaGetSymbolSize(&cluster_helper_size, cluster_real_helper));
      assert(cluster_helper_size  == sizeof(ClusterHelper));
      checkRC(cudaMemcpyToSymbol(cluster_add_matrix, &host_helper.add_matrix_, sizeof(int*))
				    );
      checkRC(cudaMemcpyToSymbol(cluster_sub_matrix, &host_helper.sub_matrix_, sizeof(int*))
				    );
      checkRC(cudaMemcpyToSymbol(cluster_real_helper, &host_helper, cluster_helper_size));
      cudaDeviceSynchronize();
    }
#ifndef NDEBUG
    checkClusterHelper<<<1, 1, 0, 0>>>(nc, lds);
#endif
  });
}

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca
