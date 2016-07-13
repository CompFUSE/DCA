// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description
//
// \Delta \tau = \frac{1}{m*N_{varpi}*step}

#ifndef MATH_LIBRARY_NFFT_DOMAINS_NFFT_TIME_DOMAIN_H
#define MATH_LIBRARY_NFFT_DOMAINS_NFFT_TIME_DOMAIN_H

#include <stdexcept>
#include <string>
#include <vector>

#include "math_library/functional_transforms/domain_specifications/domain_specifications.hpp"

namespace math_algorithms {
namespace NFFT {
// math_algorithms::NFFT::

enum NFFT_TIME_DOMAIN_NAMES { LEFT_ORIENTED, PADDED, WINDOW_FUNCTION, FOLDED_WINDOW_FUNCTION };

std::string to_str(NFFT_TIME_DOMAIN_NAMES NAME) {
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

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
class nfft_time_domain {
public:
  const static int DIMENSION = 1;

  typedef typename dnfft_type::scalar_type scalar_type;
  typedef typename dnfft_type::scalar_type element_type;

  typedef domain_specifications<scalar_type, element_type, DISCRETE, KRONECKER_DELTA, PERIODIC, EQUIDISTANT>
      discrete_periodic_dmn_1D_type;

public:
  static bool& is_initialized();

  static int& get_size();
  static std::string get_name();

  static std::vector<element_type>& get_elements();

  inline static scalar_type& first_element();

  inline static scalar_type& get_delta();
  inline static scalar_type& get_Delta();

  inline static scalar_type& get_one_div_delta();
  inline static scalar_type& get_one_div_Delta();

  static void initialize(dnfft_type& dnfft_obj);
};

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
bool& nfft_time_domain<NAME, dnfft_type>::is_initialized() {
  static bool initialized = false;
  return initialized;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
int& nfft_time_domain<NAME, dnfft_type>::get_size() {
  static int size = 0;
  return size;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
std::string nfft_time_domain<NAME, dnfft_type>::get_name() {
  std::string name = "nfft_time_domain " + to_str(NAME);
  return name;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
std::vector<typename dnfft_type::scalar_type>& nfft_time_domain<NAME, dnfft_type>::get_elements() {
  static std::vector<element_type> elements(0);
  return elements;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
typename dnfft_type::scalar_type& nfft_time_domain<NAME, dnfft_type>::first_element() {
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

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
typename dnfft_type::scalar_type& nfft_time_domain<NAME, dnfft_type>::get_delta() {
  static scalar_type delta;
  return delta;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
typename dnfft_type::scalar_type& nfft_time_domain<NAME, dnfft_type>::get_Delta() {
  static scalar_type Delta;
  return Delta;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
typename dnfft_type::scalar_type& nfft_time_domain<NAME, dnfft_type>::get_one_div_delta() {
  static scalar_type one_div_delta;
  return one_div_delta;
}

template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
typename dnfft_type::scalar_type& nfft_time_domain<NAME, dnfft_type>::get_one_div_Delta() {
  static scalar_type one_div_Delta;
  return one_div_Delta;
}

/*!
 *   \Delta = \frac{1}{OVER_SAMPLING*MAX_FREQUENCY}
 *   \delta = \frac{1}{OVER_SAMPLING*MAX_FREQUENCY}*\frac{1}{WINDOW_SAMPLING}
 */
template <NFFT_TIME_DOMAIN_NAMES NAME, typename dnfft_type>
void nfft_time_domain<NAME, dnfft_type>::initialize(dnfft_type& dnfft_obj) {
  if (not is_initialized()) {
    // cout << "\n\n\n\n\n\t\t " << __FUNCTION__ << " " << get_name() << "\n\n\n\n\n";

    is_initialized() = true;

    const int OVER_SAMPLING = dnfft_obj.get_oversampling_factor();
    const int MAX_FREQUENCY = dnfft_obj.get_maximum_frequency();
    const int WINDOW_SAMPLING = dnfft_obj.get_window_sampling_factor();

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
         *   \tau \in [-0.5-2*OVER_SAMPLING, -0.5-2*OVER_SAMPLING+delta,
         * -0.5-2*OVER_SAMPLING+2*\delta, ... , 0.5+2*OVER_SAMPLING-\delta]
         */

        get_size() = OVER_SAMPLING * MAX_FREQUENCY + 4 * OVER_SAMPLING;

        get_elements().resize(get_size(), -0.5 - 2. * OVER_SAMPLING * Delta);

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

}  // NFFT
}  // math_algorithm

#endif  // MATH_LIBRARY_NFFT_DOMAINS_NFFT_TIME_DOMAIN_H
