// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// VASP Brillouin zone cut domain.

#ifndef DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BRILLOUIN_ZONE_CUT_DOMAIN_HPP
#define DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BRILLOUIN_ZONE_CUT_DOMAIN_HPP

#include <cassert>
#include <string>
#include <vector>

#include "dca/phys/dft_connection/vasp/vasp_domains/brillouin_zone_name.hpp"
#include "dca/phys/domains/cluster/cluster_specifications.hpp"

namespace dca {
namespace phys {
namespace dft {
namespace vasp {
// dca::phys::dft::vasp::

class vasp_brillouin_zone_cut_domain {
  const static int DIMENSION = 3;

public:
  typedef domains::cluster_specifications<double, domains::VASP_LATTICE, domains::MOMENTUM_SPACE,
                                          domains::PARALLELLEPIPEDUM>::dmn_specifications_type
      dmn_specifications_type;

  typedef vasp_brillouin_zone_cut_domain this_type;
  typedef std::vector<double> element_type;

  static int& get_N() {
    static int N = 65;
    return N;
  }

  static brillouin_zone_name& get_brillouin_zone_name() {
    static brillouin_zone_name name = NO_NAME;
    return name;
  }

  static std::string& get_name() {
    static std::string name = 0;
    return name;
  }

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<std::vector<double>> v(0);
    return v;
  }

  static void initialize(std::vector<std::vector<double>>& path_vecs) {
    get_elements().resize(0);

    std::vector<double> k(DIMENSION);

    for (int i = 0; i < int(path_vecs.size()) - 1; i++) {
      std::vector<double>& p0 = path_vecs[i];
      std::vector<double>& p1 = path_vecs[i + 1];

      assert(DIMENSION == int(p0.size()));
      assert(DIMENSION == int(p1.size()));

      for (int l = 0; l < get_N(); l++) {
        for (int z = 0; z < DIMENSION; z++)
          k[z] = (1. - double(l) / double(get_N())) * p0[z] + double(l) / double(get_N()) * p1[z];

        get_elements().push_back(k);
      }
    }

    get_size() = get_elements().size();
  }
};

}  // vasp
}  // dft
}  // phys
}  // dca

#endif  // DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BRILLOUIN_ZONE_CUT_DOMAIN_HPP
