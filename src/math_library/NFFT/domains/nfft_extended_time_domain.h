//-*-C++-*-

#ifndef MATH_LIBRARY_NFFT_DOMAINS_NFFT_EXTENDED_TIME_DOMAIN_H
#define MATH_LIBRARY_NFFT_DOMAINS_NFFT_EXTENDED_TIME_DOMAIN_H

#include <vector>

namespace math_algorithms {
namespace NFFT {
// math_algorithms::NFFT::

// \Delta \tau = \frac{1}{m*N_{\varpi}}
// \tau \in [-0.5-2*m*\Delta\tau, ... , 0.5+2*m*\Delta\tau]
template <int oversampling, typename w_dmn_t>
struct nfft_extended_time_domain {
  typedef double element_type;

  static int get_size() {
    return oversampling * w_dmn_t::dmn_size() / 2 + 4 * oversampling;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements = initialize();
    return elements;
  }

private:
  static std::vector<double> initialize() {
    double delta_t = 2. / double(oversampling * w_dmn_t::dmn_size());

    std::vector<double> elements(get_size(), -0.5 - 2 * oversampling * delta_t);

    for (int l = 0; l < get_size(); l++)
      elements[l] += l * delta_t;

    return elements;
  }
};

}  // NFFT
}  // math_algorithm

#endif  // MATH_LIBRARY_NFFT_DOMAINS_NFFT_EXTENDED_TIME_DOMAIN_H
