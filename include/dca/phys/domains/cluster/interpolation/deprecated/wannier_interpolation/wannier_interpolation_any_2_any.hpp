// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class performs a Wannier interpolation on a dca::func::function.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_ANY_2_ANY_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_ANY_2_ANY_HPP

#include <cstring>

#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/interpolation/wannier_interpolation/wannier_interpolation_kernel.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename type_input, typename type_output, int dmn_number>
struct wannier_interpolation_any_2_any {
  template <typename scalartype_1, typename dmn_type_1, typename scalartype_2, typename dmn_type_2>
  static void execute(func::function<scalartype_1, dmn_type_1>& f_source,
                      func::function<scalartype_2, dmn_type_2>& f_target) {
    int Nb_sbdms = f_source.signature();
    int Nb_elements = f_source.size();

    int* coordinate = new int[Nb_sbdms];
    std::memset(coordinate, 0, sizeof(int) * Nb_sbdms);

    std::complex<double>* input_values = new std::complex<double>[f_source[dmn_number]];
    std::complex<double>* output_values = new std::complex<double>[f_target[dmn_number]];

    {
      wannier_interpolation_kernel<type_input, type_output> kernel;

      int Nb_WI = Nb_elements / f_source[dmn_number];

      for (int l = 0; l < Nb_WI; l++) {
        int linind = l;
        for (int j = Nb_sbdms - 1; j > -1; j--) {
          if (j != dmn_number) {
            coordinate[j] = linind % f_source[j];
            linind = (linind - coordinate[j]) / f_source[j];
          }
        }

        f_source.slice(dmn_number, coordinate, input_values);

        kernel.execute(input_values, output_values);

        f_target.distribute(dmn_number, coordinate, output_values);
      }
    }

    delete[] coordinate;
    delete[] input_values;
    delete[] output_values;
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_ANY_2_ANY_HPP
