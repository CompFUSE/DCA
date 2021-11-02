// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file implements cublas related utilities.

#include "dca/linalg/util/handle_functions.hpp"
#include "dca/linalg/util/util_gpublas.hpp"

#include <mutex>
#include <magma_v2.h>


namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#if defined(DCA_HAVE_GPU)
void initializeMagma() {
  static std::once_flag flag;
  std::call_once(flag, []() { magma_init(); });
}

#if defined(DCA_HAVE_CUDA)
int getGpuBLASVersion() {
  int version = 0;
  cublasStatus_t ret = cublasGetVersion(getHandle(0), &version);
  checkRC(ret);
  return version;
}
#elif defined(DCA_HAVE_HIP)
int getGpuBLASVersion() {
  int version = hipblasVersionMajor;
  return version;
} 
#endif
  
#else
void initializeMagma() {}
#endif  // DCA_HAVE_CUDA

}  // namespace util
}  // namespace linalg
}  // namespace dca
