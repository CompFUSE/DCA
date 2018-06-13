// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Folded fine time domain.

#ifndef DCA_MATH_NFFT_DOMAINS_NFFT_FOLDED_TIME_DOMAIN_HPP
#define DCA_MATH_NFFT_DOMAINS_NFFT_FOLDED_TIME_DOMAIN_HPP

#include <vector>
#include "dca/math/nfft/domains/nfft_extended_time_domain.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <int oversampling, int step, typename w_dmn_t>
struct nfft_folded_fine_time_domain {
  typedef double element_type;

  static int get_size() {
    return (4 * oversampling) * step;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements = initialize();
    return elements;
  }

private:
  static std::vector<double> initialize() {
    double t_0 = nfft_extended_time_domain<oversampling, w_dmn_t>::get_elements()[0];
    double t_1 = nfft_extended_time_domain<oversampling, w_dmn_t>::get_elements()[1];

    double delta_t = (t_1 - t_0);

    std::vector<double> elements(get_size(), -2 * oversampling * delta_t);

    int nb_col = step;
    int nb_row = 4 * oversampling;

    for (int j = 0; j < nb_col; j++)
      for (int i = 0; i < nb_row; i++)
        elements[i + nb_row * j] += double(i * step + j) * delta_t / double(step);

    return elements;
  }
};

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DOMAINS_NFFT_FOLDED_TIME_DOMAIN_HPP
