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


#if defined(DCA_HAVE_CUDA)
#include "dca/linalg/util/util_cublas.hpp"
#include "dca/linalg/util/error_cublas.hpp"
#elif defined(DCA_HAVE_HIP)
#include "dca/linalg/util/util_hipblas.hpp"
#include "dca/linalg/util/error_hipblas.hpp"
#endif

#include <mutex>
#include <magma_v2.h>
#include "dca/linalg/util/handle_functions.hpp"


namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_CUDA
int getCublasVersion() {
  int version = 0;
  cublasStatus_t ret = cublasGetVersion(getHandle(0), &version);
  checkRC(ret);
  return version;
}

void initializeMagma() {
  static std::once_flag flag;
  std::call_once(flag, []() { magma_init(); });
}
#else
void initializeMagma() {}
#endif  // DCA_HAVE_CUDA

}  // namespace util
}  // namespace linalg
}  // namespace dca
