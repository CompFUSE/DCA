//-*-C++-*-

#ifndef DCA_HYBRIDIZATION_TYPE_DEFINITIONS_H
#define DCA_HYBRIDIZATION_TYPE_DEFINITIONS_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *
     * \brief   This class defines common types for the Single-Site Hybridization Monte Carlo integrator
     * \author  Peter Staar
     * \version 1.0
     */
    template<class parameters_type, class MOMS_type>
    class MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>
    {
    public:

      /*!
       * \brief types that define the profiling
       */
      typedef typename parameters_type::concurrency_type concurrency_type;
      typedef typename parameters_type::profiler_type    profiler_type;

      /*!
       * \brief types that define the scalartype and matrix-type
       */
      typedef double                                    scalartype;
      //typedef resizeable_square_matrix<scalartype> vertex_vertex_matrix_type;
      typedef LIN_ALG::matrix<scalartype, LIN_ALG::CPU> vertex_vertex_matrix_type;

      /*!
       * \brief types that define the vertex and configuration type
       */
      typedef                              SS_CT_HYB_configuration    configuration_type;
      typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;
    };

  }

}

#endif
