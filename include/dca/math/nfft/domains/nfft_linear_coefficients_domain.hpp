// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Linear coefficients domain.

#ifndef DCA_MATH_NFFT_DOMAINS_NFFT_LINEAR_COEFFICIENTS_DOMAIN_HPP
#define DCA_MATH_NFFT_DOMAINS_NFFT_LINEAR_COEFFICIENTS_DOMAIN_HPP

#include <string>
#include <vector>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

struct nfft_linear_coefficients_domain {
  const static int DIMENSION = 1;

  typedef int scalar_type;
  typedef int element_type;

public:
  static int& get_size() {
    static int size = 2;
    return size;
  }

  static std::string get_name() {
    return "nfft_coefficients_domain";
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(0);
    return elements;
  }
};

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DOMAINS_NFFT_LINEAR_COEFFICIENTS_DOMAIN_HPP
