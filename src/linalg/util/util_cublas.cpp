// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements cublas related utilities.

#include "dca/linalg/util/util_cublas.hpp"
#include <cublas_v2.h>
#include <magma.h>
#include "dca/linalg/util/error_cublas.hpp"
#include "dca/linalg/util/handle_functions.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

int getCublasVersion() {
  int version = 0;
  cublasStatus_t ret = cublasGetVersion(getHandle(0), &version);
  checkRC(ret);
  return version;
}

void initializeMagma() {
  static bool initialized = false;
  if(!initialized) {
    magma_init();
    initialized = true;
  }
}

}  // util
}  // linalg
}  // dca
