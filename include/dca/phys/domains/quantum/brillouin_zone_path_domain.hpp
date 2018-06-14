// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the Brillouin zone path domain class that is templated on the type of the
// Brillouin zone cut.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_PATH_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_PATH_DOMAIN_HPP

#include <string>
#include <vector>

#include "dca/phys/domains/quantum/brillouin_zone_cut_type.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

//
// General class template.
//
template <BrillouinZoneCutType brillouin_zone_cut>
class brillouin_zone_path_domain {
public:
  template <class stream_type>
  static void to_JSON(stream_type& ss);

  template <class stream_type>
  static void to_JSON(stream_type& ss, std::string name, std::vector<std::vector<double>> elements);
};

//
// Full template specialization for the 2D square lattice.
//
template <>
class brillouin_zone_path_domain<SQUARE_2D_LATTICE> {
public:
  typedef std::vector<double> element_type;
  typedef brillouin_zone_path_domain<SQUARE_2D_LATTICE> this_type;

  static std::string get_name() {
    return std::string("SQUARE_2D_LATTICE");
  }

  static int get_size() {
    return get_elements().size();
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> v = initialize_elements();
    return v;
  }

private:
  const static int INTERPOLATION_NB = 64;
  static std::vector<element_type> initialize_elements();
};

//
// Full template specialization for the Fermi surface of a 2D square lattice.
//
template <>
class brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE> {
public:
  const static int INTERPOLATION_NB = 64;

  typedef std::vector<double> element_type;
  typedef brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE> this_type;

  static std::string get_name() {
    return std::string("FERMI_SURFACE_SQUARE_2D_LATTICE");
  }

  static int get_size() {
    return get_elements().size();
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> v = initialize_elements();
    return v;
  }

private:
  static std::vector<element_type> initialize_elements();
};

//
// Method definitions of the general class tempate.
// Need to be put after the declarations of the specializations since the general class template
// uses methods of the specialized classes.
//
template <BrillouinZoneCutType brillouin_zone_cut>
template <class stream_type>
void brillouin_zone_path_domain<brillouin_zone_cut>::to_JSON(stream_type& ss) {
  to_JSON(ss, brillouin_zone_path_domain<SQUARE_2D_LATTICE>::get_name(),
          brillouin_zone_path_domain<SQUARE_2D_LATTICE>::get_elements());

  ss << ",\n";
  to_JSON(ss, brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_name(),
          brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_elements());
}

template <BrillouinZoneCutType brillouin_zone_cut>
template <class stream_type>
void brillouin_zone_path_domain<brillouin_zone_cut>::to_JSON(
    stream_type& ss, std::string name, std::vector<std::vector<double>> elements) {
  ss << "\"" << name << "\" : [\n";

  for (size_t i = 0; i < elements.size(); i++) {
    ss << "[ ";
    for (size_t z = 0; z < elements[i].size(); z++) {
      if (z == elements[i].size() - 1)
        ss << elements[i][z] << "]";
      else
        ss << elements[i][z] << ", ";
    }

    if (i == elements.size() - 1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "]\n";
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_PATH_DOMAIN_HPP
