// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Doak      (doakpw@ornl.gov)
//


#ifndef DCA_LINALG_UTIL_UTIL_GPUBLAS_HPP
#define DCA_LINALG_UTIL_UTIL_GPUBLAS_HPP

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

int getGpuBLASVersion();

void initializeMagma();

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_UTIL_CUBLAS_HPP
