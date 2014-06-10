//-*-C++-*-

#ifndef FERMIONIC_ED_TYPE_DEFINITIONS_H
#define FERMIONIC_ED_TYPE_DEFINITIONS_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    struct ED_type_definitions
    {
    public:

      typedef typename parameter_type::profiler_type    profiler_t;
      typedef typename parameter_type::concurrency_type concurrency_type;

      //       typedef float       scalar_type;
      //       typedef scalar_type complex_type;

      typedef double      scalar_type;
      typedef scalar_type complex_type;

      //       typedef float                     scalar_type;
      //       typedef std::complex<scalar_type> complex_type;

      //       typedef double                    scalar_type;
      //       typedef std::complex<scalar_type> complex_type;

      typedef LIN_ALG::vector<scalar_type , LIN_ALG::CPU> vector_type;
      typedef LIN_ALG::matrix<complex_type, LIN_ALG::CPU> matrix_type;

      typedef LIN_ALG::matrix<int         , LIN_ALG::CPU> int_matrix_type;


      typedef dmn_0<occupation_number_domain<b_dmn, s_dmn, r_dmn> > occ_dmn;
      typedef dmn_0<magnetization_domain    <b_dmn, s_dmn, r_dmn> > mag_dmn;

      typedef dmn_2<b_dmn, s_dmn>     nu_dmn;
      typedef dmn_2<occ_dmn, mag_dmn> occ_mag_dmn;

      typedef dmn_3<b_dmn, s_dmn, r_dmn> b_s_r;
    };
  }

}

#endif
