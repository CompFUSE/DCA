// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements symmetrize_single_particle_function.hpp.

#include "dca/phys/dca_step/symmetrization/symmetrize_single_particle_function.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <>
void symmetrize_single_particle_function::difference(float val, std::string function_name,
                                                     std::string dmn_name) {
  if (std::abs(val) > 1.e-3) {
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t"
              << std::abs(val) << "\n\n";
    // throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <>
void symmetrize_single_particle_function::difference(float val0, float val1,
                                                     std::string function_name, std::string dmn_name) {
  if (std::abs(val0 - val1) > 1.e-3) {
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t"
              << std::abs(val0 - val1) << "\n\n";
    // throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <>
void symmetrize_single_particle_function::difference(std::complex<float> val,
                                                     std::string function_name, std::string dmn_name) {
  if (abs(val) > 1.e-3) {
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t"
              << abs(val) << "\n\n";
    // throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <>
void symmetrize_single_particle_function::difference(std::complex<float> val0,
                                                     std::complex<float> val1,
                                                     std::string function_name, std::string dmn_name) {
  if (abs(val0 - val1) > 1.e-3) {
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t"
              << abs(val0 - val1) << "\n\n";
    // throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // phys
}  // dca
