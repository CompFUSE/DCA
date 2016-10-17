// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the NFFT time domain names.

#ifndef DCA_MATH_NFFT_DOMAINS_NFFT_TIME_DOMAIN_NAMES_HPP
#define DCA_MATH_NFFT_DOMAINS_NFFT_TIME_DOMAIN_NAMES_HPP

#include <string>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

enum NfftTimeDomainNames { LEFT_ORIENTED, PADDED, WINDOW_FUNCTION, FOLDED_WINDOW_FUNCTION };

std::string to_str(NfftTimeDomainNames NAME);

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DOMAINS_NFFT_TIME_DOMAIN_HPP
