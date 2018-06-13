// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// \Delta \tau = \frac{1}{m*N_{varpi}*step}

#ifndef DCA_MATH_NFFT_DOMAINS_NFFT_TIME_DOMAIN_HPP
#define DCA_MATH_NFFT_DOMAINS_NFFT_TIME_DOMAIN_HPP

#include <stdexcept>
#include <string>
#include <vector>

#include "dca/math/function_transform/domain_specifications.hpp"
#include "dca/math/nfft/domains/nfft_time_domain_names.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <NfftTimeDomainNames NAME, typename dnfft_type>
class nfft_time_domain {
public:
  const static int DIMENSION = 1;

  using scalar_type = typename dnfft_type::ElementType;
  using element_type = typename dnfft_type::ElementType;

  typedef transform::domain_specifications<scalar_type, element_type, transform::DISCRETE,
                                           transform::KRONECKER_DELTA, transform::PERIODIC,
                                           transform::EQUIDISTANT>
      discrete_periodic_dmn_1D_type;

  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    std::string name = "nfft_time_domain " + to_str(NAME);
    return name;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(0);
    return elements;
  }

  static scalar_type& first_element();

  static scalar_type& get_delta() {
    static scalar_type delta;
    return delta;
  }
  static scalar_type& get_Delta() {
    static scalar_type Delta;
    return Delta;
  }

  static scalar_type& get_one_div_delta() {
    static scalar_type one_div_delta;
    return one_div_delta;
  }

  static scalar_type& get_one_div_Delta() {
    static scalar_type one_div_Delta;
    return one_div_Delta;
  }

  static void initialize(const dnfft_type& dnfft_obj);
};

template <NfftTimeDomainNames NAME, typename dnfft_type>
typename nfft_time_domain<NAME, dnfft_type>::scalar_type& nfft_time_domain<NAME, dnfft_type>::first_element() {
  switch (NAME) {
    case LEFT_ORIENTED: {
      static scalar_type first;
      return first;
    } break;

    case PADDED: {
      static scalar_type first;
      return first;
    } break;

    case WINDOW_FUNCTION: {
      static scalar_type first;
      return first;
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

/*!
 *   \Delta = \frac{1}{OVER_SAMPLING*MAX_FREQUENCY}
 *   \delta = \frac{1}{OVER_SAMPLING*MAX_FREQUENCY}*\frac{1}{WINDOW_SAMPLING}
 */
template <NfftTimeDomainNames NAME, typename dnfft_type>
void nfft_time_domain<NAME, dnfft_type>::initialize(const dnfft_type& dnfft_obj) {
  if (not is_initialized()) {
    // cout << "\n\n\n\n\n\t\t " << __FUNCTION__ << " " << get_name() << "\n\n\n\n\n";

    is_initialized() = true;

    const int OVER_SAMPLING = dnfft_obj.get_oversampling();
    const int MAX_FREQUENCY = dnfft_obj.maximumFrequency();
    const int WINDOW_SAMPLING = dnfft_obj.get_window_sampling();

    scalar_type Delta = 1. / scalar_type(OVER_SAMPLING * MAX_FREQUENCY);
    scalar_type delta = 1. / scalar_type(OVER_SAMPLING * MAX_FREQUENCY * WINDOW_SAMPLING);

    get_delta() = delta;
    get_Delta() = Delta;

    get_one_div_delta() = 1. / delta;
    get_one_div_Delta() = 1. / Delta;

    switch (NAME) {
      case LEFT_ORIENTED: {
        /*!
         *   \tau \in [-0.5, -0.5+delta, -0.5+2*\delta, ... , 0.5-\delta]
         */

        get_size() = OVER_SAMPLING * MAX_FREQUENCY;

        get_elements().resize(get_size(), -0.5);

        for (int l = 0; l < get_size(); l++)
          get_elements()[l] += l * Delta;

        first_element() = get_elements()[0];
      } break;

      case PADDED: {
        /*!
         *   \tau \in [-0.5-OVER_SAMPLING*\delta, -0.5-OVER_SAMPLING*\delta+\delta,
         *             -0.5-OVER_SAMPLING*\delta+2*\delta, ... , 0.5+OVER_SAMPLING*\delta-\delta]
         */

        get_size() = OVER_SAMPLING * MAX_FREQUENCY + 2 * OVER_SAMPLING;

        get_elements().resize(get_size(), -0.5 - OVER_SAMPLING * Delta);

        for (int l = 0; l < get_size(); l++)
          get_elements()[l] += l * Delta;

        first_element() = get_elements()[0];
      } break;

      case WINDOW_FUNCTION: {
        /*!
         *   \tau \in [-2./scalar_type(MAX_FREQUENCY), ... , 2./scalar_type(MAX_FREQUENCY)]
         */

        get_size() = 4 * OVER_SAMPLING * WINDOW_SAMPLING;

        get_elements().resize(get_size(), -2. / scalar_type(MAX_FREQUENCY));

        for (int l = 0; l < get_size(); l++)
          get_elements()[l] += l * delta;

        first_element() = get_elements()[0];
      } break;

      case FOLDED_WINDOW_FUNCTION: {
        /*!
         *   \tau \in [-0.5-2*OVER_SAMPLING, -0.5-2*OVER_SAMPLING+delta,
         * -0.5-2*OVER_SAMPLING+2*\delta, ... , 0.5+2*OVER_SAMPLING-\delta]
         */

        get_size() = 4 * OVER_SAMPLING * WINDOW_SAMPLING;

        get_elements().resize(get_size(), -2. / scalar_type(MAX_FREQUENCY));

        int nb_row = 4 * OVER_SAMPLING;
        int nb_col = WINDOW_SAMPLING;

        for (int j = 0; j < nb_col; j++)
          for (int i = 0; i < nb_row; i++)
            get_elements()[i + nb_row * j] += scalar_type(i * WINDOW_SAMPLING + j) * delta;
      } break;

      default:
        throw std::logic_error(__FUNCTION__);
    }
  }
}

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DOMAINS_NFFT_TIME_DOMAIN_HPP
