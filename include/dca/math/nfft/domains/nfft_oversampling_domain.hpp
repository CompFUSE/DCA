// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Oversampling domain.

#ifndef DCA_MATH_NFFT_DOMAINS_NFFT_OVERSAMPLING_DOMAIN_HPP
#define DCA_MATH_NFFT_DOMAINS_NFFT_OVERSAMPLING_DOMAIN_HPP

#include <string>
#include <vector>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <typename dnfft_type>
struct nfft_oversampling_domain {
  const static int DIMENSION = 1;

  typedef int scalar_type;
  typedef int element_type;

public:
  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    return "nfft_oversampling_domain";
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(0);
    return elements;
  }

  static void initialize(const dnfft_type& dnfft_obj) {
    get_size() = 4 * dnfft_obj.get_oversampling();
  }
};

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DOMAINS_NFFT_OVERSAMPLING_DOMAIN_HPP
