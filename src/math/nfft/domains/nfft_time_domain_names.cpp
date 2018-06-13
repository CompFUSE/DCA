// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements nfft_time_domain_names.cpp.

#include "dca/math/nfft/domains/nfft_time_domain_names.hpp"
#include <stdexcept>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

std::string to_str(NfftTimeDomainNames NAME) {
  switch (NAME) {
    case LEFT_ORIENTED:
      return "LEFT_ORIENTED";
      break;

    case PADDED:
      return "PADDED";
      break;

    case WINDOW_FUNCTION:
      return "WINDOW_FUNCTION";
      break;

    case FOLDED_WINDOW_FUNCTION:
      return "FOLDED_WINDOW_FUNCTION";
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // nfft
}  // math
}  // dca
