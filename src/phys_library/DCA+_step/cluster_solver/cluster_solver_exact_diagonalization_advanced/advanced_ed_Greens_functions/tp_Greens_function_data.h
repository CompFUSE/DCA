//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_TP_GREENS_FUNCTION_DATA_H
#define ADVANCED_FERMIONIC_TP_GREENS_FUNCTION_DATA_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {
    template<typename ed_options>
    class tp_Greens_function_data
    {

      typedef typename ed_options::b_dmn b_dmn;
      typedef typename ed_options::s_dmn s_dmn;
      typedef typename ed_options::r_dmn r_dmn;
      typedef typename ed_options::k_dmn k_dmn;

      typedef typename ed_options::profiler_t       profiler_t;
      typedef typename ed_options::concurrency_type concurrency_type;

      typedef typename ed_options::scalar_type  scalar_type;
      typedef typename ed_options::complex_type complex_type;

      typedef typename ed_options::vector_type         vector_type;
      typedef typename ed_options::matrix_type         matrix_type;
      typedef typename ed_options::int_matrix_type int_matrix_type;

      typedef typename ed_options::nu_dmn        nu_dmn;
      typedef typename ed_options::nu_r_dmn_type nu_r_dmn_type;

      typedef typename ed_options::b_s_r  b_s_r_dmn_type;

      typedef typename ed_options::bs_dmn_type  bs_dmn_type;
      typedef typename ed_options::bsr_dmn_type bsr_dmn_type;

      typedef typename ed_options::nu_nu_r_dmn_type nu_nu_r_dmn_type;

      typedef tp_Greens_function_data<ed_options> this_type;

      typedef dmn_2<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED> wm_wn_dmn_type;

      typedef dmn_7<nu, nu, nu, nu, r_dmn, r_dmn, r_dmn> nu_nu_nu_nu_r_r_r_dmn_type;      

    public:

      tp_Greens_function_data();

      tp_Greens_function_data(const this_type& other);

      ~tp_Greens_function_data();

      void set_indices(int l);

      template<typename parameter_type>
      void initialize(parameter_type& parameters);

      void sum_to(FUNC_LIB::function<complex_type, dmn_4<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, nu_nu_nu_nu_r_r_r_dmn_type> >& G_tp_ref);

    public:

      int b_i;
      int s_i;

      int b_j;
      int s_j;

      int r_i;
      int r_j;

      int delta_r;

      int nu_i;
      int nu_j;

      int bsr_i;
      int bsr_j;

      bs_dmn_type  bs_dmn;
      bsr_dmn_type bsr_dmn;

      nu_r_dmn_type         nu_r_dmn;
      nu_nu_r_dmn_type nu_nu_r_dmn;

      nu_nu_nu_nu_r_r_r_dmn_type nu_nu_nu_nu_r_r_r_dmn;

      matrix_type tmp;

      matrix_type overlap_0;
      matrix_type overlap_1;
      matrix_type overlap_2;
      matrix_type overlap_3;

      FUNC_LIB::function<complex_type, dmn_4<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, nu_nu_nu_nu_r_r_r_dmn_type> > G_tp;
    };

    template<typename ed_options>
    tp_Greens_function_data<ed_options>::tp_Greens_function_data()
    {}

    template<typename ed_options>
    tp_Greens_function_data<ed_options>::tp_Greens_function_data(const this_type& /*other*/)
    {}
    
    template<typename ed_options>
    tp_Greens_function_data<ed_options>::~tp_Greens_function_data()
    {}

    template<typename ed_options>
    template<typename parameter_type>
    void tp_Greens_function_data<ed_options>::initialize(parameter_type& /*parameters*/)
    {
      G_tp = 0;
    }

    template<typename ed_options>
    void tp_Greens_function_data<ed_options>::sum_to(FUNC_LIB::function<complex_type, dmn_4<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, nu_nu_nu_nu_r_r_r_dmn_type> >& G_tp_ref)
    {
      G_tp_ref += G_tp;
    }

  }

}

#endif
