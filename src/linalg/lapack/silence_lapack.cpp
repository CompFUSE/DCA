// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements silence_lapack.hpp.

#include "dca/linalg/lapack/silence_lapack.hpp"

#ifdef DCA_WITH_ESSL
#include <essl.h>
#endif  // DCA_WITH_ESSL

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

#ifdef DCA_WITH_ESSL
void silenceLapack() {
  const _ESVINT code = 2146;
  static _ESVINT info1(0), info2(0);
  einfo(code, info1, info2);

  static _ESVINT iusadr(1);
  errset(code, 256, -1, 0, static_cast<const void*>(&iusadr), 0);
}

#else   // DCA_WITH_ESSL
void silenceLapack() {}
#endif  // DCA_WITH_ESSL

}  // lapack
}  // linalg
}  // dca
