//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_ED_TYPE_DEFINITIONS_H
#define ADVANCED_FERMIONIC_ED_TYPE_DEFINITIONS_H

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {

    enum phi_names {PHI_SINGLET, PHI_MULTIPLET};

    //typedef phi_names 
    //enum psi_names {PSI_SINGLET, PSI_MULTIPLET};

    template<typename parameter_type>
    struct advanced_ed_options
    {
#include "type_definitions.h"

    public:

      const static size_t N=8*sizeof(size_t);

      typedef std::bitset<N> phi_type;

      typedef typename parameter_type::profiler_type    profiler_t;
      typedef typename parameter_type::concurrency_type concurrency_type;

      typedef int                       int_type;

      //typedef float                    scalar_type;
      typedef double                    scalar_type;
      typedef std::complex<scalar_type> complex_type;

      typedef LIN_ALG::vector<scalar_type , LIN_ALG::CPU> vector_type;
      typedef LIN_ALG::matrix<complex_type, LIN_ALG::CPU> matrix_type;

      typedef LIN_ALG::matrix<int         , LIN_ALG::CPU> int_matrix_type;

      typedef b     b_dmn;
      typedef s     s_dmn;
      typedef r_DCA r_dmn;
      typedef k_DCA k_dmn;

      typedef dmn_2<b_dmn, s_dmn>        nu_dmn;
      typedef dmn_2<b_dmn, s_dmn>        bs_dmn_type;

      typedef dmn_2<nu_dmn, r_dmn>       nu_r_dmn_type;

      typedef dmn_3<b_dmn, s_dmn, r_dmn> b_s_r;
      typedef dmn_3<b_dmn, s_dmn, r_dmn> bsr_dmn_type;

      typedef dmn_3<b_dmn, s_dmn, k_dmn> b_s_k;
      typedef dmn_3<b_dmn, s_dmn, k_dmn> bsk_dmn_type;

      typedef dmn_3<nu_dmn, nu_dmn, r_dmn> nu_nu_r_dmn_type;
      typedef dmn_3<nu_dmn, nu_dmn, k_dmn> nu_nu_k_dmn_type;

    public:

      static scalar_type get_epsilon()
      {
        return 1.e-3;
      }

    };

  }

}

#endif
