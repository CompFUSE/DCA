// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the nonlocal \f$\chi(k_1,k_2,q)\f$ in the bands.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_NONLOCAL_CHI_ATOMIC_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_NONLOCAL_CHI_ATOMIC_HPP

#include <complex>
#include <stdexcept>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

//
// Empty class template.
//
template <class model_type, FourPointType four_point>
class accumulator_nonlocal_chi_atomic {};

//
// Specialization for particle-hole-magnetic channel.
//
template <class model_type>
class accumulator_nonlocal_chi_atomic<model_type, PARTICLE_HOLE_MAGNETIC> {
public:
  using b = func::dmn_0<domains::electron_band_domain>;

  template <typename scalartype>
  static inline void execute(std::complex<scalartype>* G2_DN_n1_m2_k1_k2_w1_w2,
                             std::complex<scalartype>* G2_UP_n1_m2_k1_k2_w1_w2,
                             std::complex<scalartype>* G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                             std::complex<scalartype>* G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                             std::complex<scalartype>* G2_DN_n1_m1_k1_k1_plus_q_w1_w1,
                             std::complex<scalartype>* G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
                             std::complex<scalartype>* G2_DN_n2_m2_k2_plus_q_k2_w2_w2,
                             std::complex<scalartype>* G2_UP_n2_m2_k2_plus_q_k2_w2_w2,
                             std::complex<double>* G4_k_k_w_w, scalartype& sign) {
    const static int BANDS = b::dmn_size();

    std::complex<scalartype> G4;

    for (int m1 = 0; m1 < BANDS; m1++) {
      for (int m2 = 0; m2 < BANDS; m2++) {
        for (int n1 = 0; n1 < BANDS; n1++) {
          for (int n2 = 0; n2 < BANDS; n2++) {
            G4 = -(G2_DN_n1_m2_k1_k2_w1_w2[n1 + BANDS * m2] *
                       G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2 + BANDS * m1] +
                   G2_UP_n1_m2_k1_k2_w1_w2[n1 + BANDS * m2] *
                       G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2 + BANDS * m1])

                 +
                 (G2_UP_n1_m1_k1_k1_plus_q_w1_w1[n1 + BANDS * m1] -
                  G2_DN_n1_m1_k1_k1_plus_q_w1_w1[n1 + BANDS * m1]) *
                     (G2_UP_n2_m2_k2_plus_q_k2_w2_w2[n2 + BANDS * m2] -
                      G2_DN_n2_m2_k2_plus_q_k2_w2_w2[n2 + BANDS * m2]);

            /*
              G4 = - (G2_k_k_w_w_e_DN(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_DN(n2, m1, k2_plus_q,
              k1_plus_q, w2+w_channel, w1+w_channel)
              + G2_k_k_w_w_e_UP(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_UP(n2, m1, k2_plus_q,
              k1_plus_q, w2+w_channel, w1+w_channel))

              + (G2_k_k_w_w_e_UP(n1, m1, k1, k1_plus_q, w1, w1+w_channel) - G2_k_k_w_w_e_DN(n1, m1,
              k1, k1_plus_q, w1, w1+w_channel))
              * (G2_k_k_w_w_e_UP(n2, m2, k2_plus_q, k2, w2+w_channel, w2) - G2_k_k_w_w_e_DN(n2, m2,
              k2_plus_q, k2, w2+w_channel, w2));
              */

            G4_k_k_w_w[n1 + BANDS * n2 + BANDS * BANDS * m1 + BANDS * BANDS * BANDS * m2] +=
                std::complex<double>(sign * G4);
          }
        }
      }
    }
  }
};

//
// Specialization for particle-hole-magnetic channel and tight-binding model.
//
template <class lattice_t>
class accumulator_nonlocal_chi_atomic<dca::phys::models::TightBindingModel<lattice_t>,
                                      PARTICLE_HOLE_MAGNETIC> {
  const static int BANDS = lattice_t::BANDS;

public:
  template <typename scalartype>
  static inline void execute(std::complex<scalartype>* G2_DN_n1_m2_k1_k2_w1_w2,
                             std::complex<scalartype>* G2_UP_n1_m2_k1_k2_w1_w2,
                             std::complex<scalartype>* G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                             std::complex<scalartype>* G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                             std::complex<scalartype>* G2_DN_n1_m1_k1_k1_plus_q_w1_w1,
                             std::complex<scalartype>* G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
                             std::complex<scalartype>* G2_DN_n2_m2_k2_plus_q_k2_w2_w2,
                             std::complex<scalartype>* G2_UP_n2_m2_k2_plus_q_k2_w2_w2,
                             std::complex<double>* G4_k_k_w_w, scalartype& sign) {
    std::complex<scalartype> G4;

    for (int m1 = 0; m1 < BANDS; m1++) {
      for (int m2 = 0; m2 < BANDS; m2++) {
        for (int n1 = 0; n1 < BANDS; n1++) {
          for (int n2 = 0; n2 < BANDS; n2++) {
            G4 = -(G2_DN_n1_m2_k1_k2_w1_w2[n1 + BANDS * m2] *
                       G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2 + BANDS * m1] +
                   G2_UP_n1_m2_k1_k2_w1_w2[n1 + BANDS * m2] *
                       G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2 + BANDS * m1])

                 +
                 (G2_UP_n1_m1_k1_k1_plus_q_w1_w1[n1 + BANDS * m1] -
                  G2_DN_n1_m1_k1_k1_plus_q_w1_w1[n1 + BANDS * m1]) *
                     (G2_UP_n2_m2_k2_plus_q_k2_w2_w2[n2 + BANDS * m2] -
                      G2_DN_n2_m2_k2_plus_q_k2_w2_w2[n2 + BANDS * m2]);

            /*
              G4 = - (G2_k_k_w_w_e_DN(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_DN(n2, m1, k2_plus_q,
              k1_plus_q, w2+w_channel, w1+w_channel)
              + G2_k_k_w_w_e_UP(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_UP(n2, m1, k2_plus_q,
              k1_plus_q, w2+w_channel, w1+w_channel))

              + (G2_k_k_w_w_e_UP(n1, m1, k1, k1_plus_q, w1, w1+w_channel) - G2_k_k_w_w_e_DN(n1, m1,
              k1, k1_plus_q, w1, w1+w_channel))
              * (G2_k_k_w_w_e_UP(n2, m2, k2_plus_q, k2, w2+w_channel, w2) - G2_k_k_w_w_e_DN(n2, m2,
              k2_plus_q, k2, w2+w_channel, w2));
              */

            G4_k_k_w_w[n1 + BANDS * n2 + BANDS * BANDS * m1 + BANDS * BANDS * BANDS * m2] +=
                std::complex<double>(sign * G4);
          }
        }
      }
    }
  }
};

//
// Specialization for particle-hole-charge channel.
//
template <class model_type>
class accumulator_nonlocal_chi_atomic<model_type, PARTICLE_HOLE_CHARGE> {
  inline void execute() {
    throw std::logic_error(__FUNCTION__);
  }
};

//
// Specialization for particle-hole-superconducting channel.
//
template <class model_type>
class accumulator_nonlocal_chi_atomic<model_type, PARTICLE_PARTICLE_UP_DOWN> {
  inline void execute() {
    throw std::logic_error(__FUNCTION__);
  }
};

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_NONLOCAL_CHI_ATOMIC_HPP
