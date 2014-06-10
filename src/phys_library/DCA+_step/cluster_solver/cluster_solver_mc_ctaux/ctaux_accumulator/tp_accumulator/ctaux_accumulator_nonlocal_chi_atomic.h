//-*-C++-*-

#ifndef DCA_QMCI_ACCUMULATOR_NONLOCAL_CHI_ATOMIC_H
#define DCA_QMCI_ACCUMULATOR_NONLOCAL_CHI_ATOMIC_H

namespace DCA
{
  namespace QMCI
  {
    namespace CT_AUX_ACCUMULATION
    {
      /*!
       *  \class   accumulator_nonlocal_chi_atomic
       *  \ingroup CT-AUX-ACCUMULATOR
       *
       *  \brief   this class computes the nonlocal \f$\chi(k_1,k_2,q)\f$ in the bands
       *  \date    12-07-2011
       *  \author  Peter Staar
       *  \version 1.0
       */
      template<class model_type, vertex_measurement_type vertex_measurement>
      class accumulator_nonlocal_chi_atomic
      {};

      template<class model_type>
      class accumulator_nonlocal_chi_atomic<model_type, PARTICLE_HOLE_MAGNETIC>
      {
#include "type_definitions.h"

      public:

        template<typename scalartype>
        static inline void execute(std::complex<scalartype>* G2_DN_n1_m2_k1_k2_w1_w2              , std::complex<scalartype>* G2_UP_n1_m2_k1_k2_w1_w2,
                                   std::complex<scalartype>* G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, std::complex<scalartype>* G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                                   std::complex<scalartype>* G2_DN_n1_m1_k1_k1_plus_q_w1_w1       , std::complex<scalartype>* G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
                                   std::complex<scalartype>* G2_DN_n2_m2_k2_plus_q_k2_w2_w2       , std::complex<scalartype>* G2_UP_n2_m2_k2_plus_q_k2_w2_w2,
                                   std::complex<double>* G4_k_k_w_w, scalartype& sign)
        {
          const static int BANDS = b::dmn_size();

          std::complex<scalartype> G4;

          for(int m1=0; m1<BANDS; m1++){
            for(int m2=0; m2<BANDS; m2++){
              for(int n1=0; n1<BANDS; n1++){
                for(int n2=0; n2<BANDS; n2++){

                  G4 = - (G2_DN_n1_m2_k1_k2_w1_w2[n1+BANDS*m2] * G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2+BANDS*m1]
                          + G2_UP_n1_m2_k1_k2_w1_w2[n1+BANDS*m2] * G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2+BANDS*m1])

                    + (G2_UP_n1_m1_k1_k1_plus_q_w1_w1[n1+BANDS*m1] - G2_DN_n1_m1_k1_k1_plus_q_w1_w1[n1+BANDS*m1])
                    * (G2_UP_n2_m2_k2_plus_q_k2_w2_w2[n2+BANDS*m2] - G2_DN_n2_m2_k2_plus_q_k2_w2_w2[n2+BANDS*m2]);

                  /*
                    G4 = - (G2_k_k_w_w_e_DN(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_DN(n2, m1, k2_plus_q, k1_plus_q, w2+w_channel, w1+w_channel)
                    + G2_k_k_w_w_e_UP(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_UP(n2, m1, k2_plus_q, k1_plus_q, w2+w_channel, w1+w_channel))

                    + (G2_k_k_w_w_e_UP(n1, m1, k1, k1_plus_q, w1, w1+w_channel) - G2_k_k_w_w_e_DN(n1, m1, k1, k1_plus_q, w1, w1+w_channel))
                    * (G2_k_k_w_w_e_UP(n2, m2, k2_plus_q, k2, w2+w_channel, w2) - G2_k_k_w_w_e_DN(n2, m2, k2_plus_q, k2, w2+w_channel, w2));
                    */

                  G4_k_k_w_w[n1+BANDS*n2+BANDS*BANDS*m1+BANDS*BANDS*BANDS*m2] += std::complex<double>(sign * G4);
                }
              }
            }
          }
        }
      };

      template<class lattice_t, class interaction_t>
      class accumulator_nonlocal_chi_atomic<tight_binding_model<lattice_t, interaction_t>, PARTICLE_HOLE_MAGNETIC>
      {
#include "type_definitions.h"

        const static int BANDS = lattice_t::BANDS;

      public:

        template<typename scalartype>
        static inline void execute(std::complex<scalartype>* G2_DN_n1_m2_k1_k2_w1_w2              , std::complex<scalartype>* G2_UP_n1_m2_k1_k2_w1_w2,
                                   std::complex<scalartype>* G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, std::complex<scalartype>* G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                                   std::complex<scalartype>* G2_DN_n1_m1_k1_k1_plus_q_w1_w1       , std::complex<scalartype>* G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
                                   std::complex<scalartype>* G2_DN_n2_m2_k2_plus_q_k2_w2_w2       , std::complex<scalartype>* G2_UP_n2_m2_k2_plus_q_k2_w2_w2,
                                   std::complex<double>* G4_k_k_w_w, scalartype& sign)
        {
          std::complex<scalartype> G4;

          for(int m1=0; m1<BANDS; m1++){
            for(int m2=0; m2<BANDS; m2++){
              for(int n1=0; n1<BANDS; n1++){
                for(int n2=0; n2<BANDS; n2++){

                  G4 = - (G2_DN_n1_m2_k1_k2_w1_w2[n1+BANDS*m2] * G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2+BANDS*m1]
                          + G2_UP_n1_m2_k1_k2_w1_w2[n1+BANDS*m2] * G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1[n2+BANDS*m1])

                    + (G2_UP_n1_m1_k1_k1_plus_q_w1_w1[n1+BANDS*m1] - G2_DN_n1_m1_k1_k1_plus_q_w1_w1[n1+BANDS*m1])
                    * (G2_UP_n2_m2_k2_plus_q_k2_w2_w2[n2+BANDS*m2] - G2_DN_n2_m2_k2_plus_q_k2_w2_w2[n2+BANDS*m2]);

                  /*
                    G4 = - (G2_k_k_w_w_e_DN(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_DN(n2, m1, k2_plus_q, k1_plus_q, w2+w_channel, w1+w_channel)
                    + G2_k_k_w_w_e_UP(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_UP(n2, m1, k2_plus_q, k1_plus_q, w2+w_channel, w1+w_channel))

                    + (G2_k_k_w_w_e_UP(n1, m1, k1, k1_plus_q, w1, w1+w_channel) - G2_k_k_w_w_e_DN(n1, m1, k1, k1_plus_q, w1, w1+w_channel))
                    * (G2_k_k_w_w_e_UP(n2, m2, k2_plus_q, k2, w2+w_channel, w2) - G2_k_k_w_w_e_DN(n2, m2, k2_plus_q, k2, w2+w_channel, w2));
                    */

                  G4_k_k_w_w[n1+BANDS*n2+BANDS*BANDS*m1+BANDS*BANDS*BANDS*m2] += std::complex<double>(sign * G4);
                }
              }
            }
          }
        }
      };

      template<class model_type>
      class accumulator_nonlocal_chi_atomic<model_type, PARTICLE_HOLE_CHARGE>
      {
        inline void execute()
        {
          throw std::logic_error(__FUNCTION__);
        }
      };


      template<class model_type>
      class accumulator_nonlocal_chi_atomic<model_type, PARTICLE_PARTICLE_SUPERCONDUCTING>
      {
        inline void execute()
        {
          throw std::logic_error(__FUNCTION__);
        }
      };


    }

  }

}

#endif

