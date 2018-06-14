// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Window sampling domain.

#ifndef DCA_MATH_NFFT_DOMAINS_NFFT_WINDOW_SAMPLING_DOMAIN_HPP
#define DCA_MATH_NFFT_DOMAINS_NFFT_WINDOW_SAMPLING_DOMAIN_HPP

#include <string>
#include <vector>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <typename dnfft_type>
struct nfft_window_sampling_domain {
  const static int DIMENSION = 1;

  typedef int scalar_type;
  typedef int element_type;

public:
  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    return "nfft_window_sampling_domain";
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(0);
    return elements;
  }

  static void initialize(const dnfft_type& dnfft_obj) {
    get_size() = dnfft_obj.get_window_sampling();
  }
};

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DOMAINS_NFFT_WINDOW_SAMPLING_DOMAIN_HPP
