//-*-C++-*-
// Author: Peter Staar
#ifndef MATH_LIBRARY_DOMAIN_SPECIFICATIONS_GENERIC_DOMAIN_H
#define MATH_LIBRARY_DOMAIN_SPECIFICATIONS_GENERIC_DOMAIN_H

#include <string>
#include <vector>

namespace math_algorithms {

enum DOMAIN_NAMES { DMN_0, DMN_1, DMN_2, DMN_3 };

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
class generic_domain {
public:
  typedef dmn_specs_type dmn_specifications_type;

  typedef typename dmn_specifications_type::scalar_type scalar_type;
  typedef typename dmn_specifications_type::element_type element_type;

public:
  static bool& is_initialized();

  static int& get_size();

  static int*& get_dimensions();

  static scalar_type*& get_min();
  static scalar_type*& get_max();

  static std::string& get_name();

  static std::vector<element_type>& get_elements();

  static void reset();
};

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
bool& generic_domain<NAME, dmn_specs_type>::is_initialized() {
  static bool initialized = false;
  return initialized;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
int& generic_domain<NAME, dmn_specs_type>::get_size() {
  static int size = 0;
  return size;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
int*& generic_domain<NAME, dmn_specs_type>::get_dimensions() {
  static int* dimensions = NULL;
  return dimensions;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
typename generic_domain<NAME, dmn_specs_type>::scalar_type*& generic_domain<NAME,
                                                                            dmn_specs_type>::get_min() {
  static scalar_type* min = NULL;
  return min;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
typename generic_domain<NAME, dmn_specs_type>::scalar_type*& generic_domain<NAME,
                                                                            dmn_specs_type>::get_max() {
  static scalar_type* max = NULL;
  return max;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
std::string& generic_domain<NAME, dmn_specs_type>::get_name() {
  static std::string name = dmn_specifications_type::to_str();
  return name;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
std::vector<typename generic_domain<NAME, dmn_specs_type>::element_type>& generic_domain<
    NAME, dmn_specs_type>::get_elements() {
  static std::vector<element_type> elements(get_size());
  return elements;
}

template <DOMAIN_NAMES NAME, typename dmn_specs_type>
void generic_domain<NAME, dmn_specs_type>::reset() {
  get_size() = 0;
  get_name() = "";

  get_elements().resize(0);

  is_initialized() = false;
}

}  // math_algorithms

#endif  // MATH_LIBRARY_DOMAIN_SPECIFICATIONS_GENERIC_DOMAIN_H
