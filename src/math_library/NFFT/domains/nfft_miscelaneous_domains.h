//-*-C++-*-

#ifndef MATH_LIBRARY_NFFT_DOMAINS_NFFT_MISCELANEOUS_DOMAINS_H
#define MATH_LIBRARY_NFFT_DOMAINS_NFFT_MISCELANEOUS_DOMAINS_H

#include <vector>

namespace math_algorithms {
namespace NFFT {

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

struct nfft_cubic_coefficients_domain {
  const static int DIMENSION = 1;

  typedef int scalar_type;
  typedef int element_type;

public:
  static int& get_size() {
    static int size = 4;
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

  static void initialize(dnfft_type& dnfft_obj) {
    get_size() = 4 * dnfft_obj.get_oversampling_factor();
  }
};

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

  static void initialize(dnfft_type& dnfft_obj) {
    get_size() = dnfft_obj.get_window_sampling_factor();
  }
};

}  // NFFT
}  // math_algorithm

#endif  // MATH_LIBRARY_NFFT_DOMAINS_NFFT_MISCELANEOUS_DOMAINS_H
