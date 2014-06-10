//-*-C++-*-

#ifndef COMPUTE_FIRST_ORDER_SIGMA_H
#define COMPUTE_FIRST_ORDER_SIGMA_H

namespace DCA
{

  namespace SERIES_EXPANSION
  {
    /*!
     * \author Peter Staar
     *
     * \brief  This class implements the computation of the self-energy in second order.
     */
    template<class parameter_type, class k_dmn_t>
    class sigma_perturbation<1, parameter_type, k_dmn_t>
    {
#include "type_definitions.h"

    public:

      typedef typename compute_interaction::function_type U_type;

      typedef function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > function_type;

    public:

      sigma_perturbation(parameter_type&      parameters_ref,
                         compute_interaction& interaction_obj);

      ~sigma_perturbation();

      function_type& get_function() { return Sigma; }

      void execute_on_cluster(function<double, nu>& occupancy);

      template<IO::FORMAT DATA_FORMAT>
      void write(IO::writer<DATA_FORMAT>& writer);

    protected:

      parameter_type& parameters;

      U_type& U;

      function<std::complex<double>, dmn_2<nu,nu> > d_matrix;
      function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > Sigma;
    };

    template<class parameter_type, class k_dmn_t>
    sigma_perturbation<1, parameter_type, k_dmn_t>::sigma_perturbation(parameter_type&      parameters_ref,
                                                                       compute_interaction& interaction_obj):
      parameters(parameters_ref),

      U(interaction_obj.get_function()),

      d_matrix("d-matrix"),
      Sigma("Sigma-1st-order")
    {}

    template<class parameter_type, class k_dmn_t>
    sigma_perturbation<1, parameter_type, k_dmn_t>::~sigma_perturbation()
    {}

    template<class parameter_type, class k_dmn_t>
    template<IO::FORMAT DATA_FORMAT>
    void sigma_perturbation<1, parameter_type, k_dmn_t>::write(IO::writer<DATA_FORMAT>& writer)
    {

    }

    template<class parameter_type, class k_dmn_t>
    void sigma_perturbation<1, parameter_type, k_dmn_t>::execute_on_cluster(function<double, nu>& occupancy)
    {
      Sigma = 0.;

      for(int w_ind=0; w_ind<w::dmn_size(); ++w_ind){
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

          for(int i=0; i<b::dmn_size(); ++i)
            for(int si=0; si<s::dmn_size(); ++si)
              for(int j=0; j<b::dmn_size(); ++j)
                for(int sj=0; sj<s::dmn_size(); ++sj)
                  Sigma(i, si, i, si, k_ind, w_ind) += U(i,si, j,sj)*(occupancy(j,sj)-1./2.);
        }
      }
    }


  }

}

#endif
