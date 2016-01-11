//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_SP_GREENS_FUNCTION_DATA_H
#define ADVANCED_FERMIONIC_SP_GREENS_FUNCTION_DATA_H

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {
    struct c_operator
    {
      int  index;
      int  bsr_ind;
      bool creation;
    };

    template<typename ed_options>
    class sp_Greens_function_data
    {
#include "type_definitions.h"

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

      typedef sp_Greens_function_data<ed_options> this_type;

    public:

      sp_Greens_function_data();

      sp_Greens_function_data(const this_type& other);

      ~sp_Greens_function_data();

      //       this_type& operator=(      this_type& other);
      //       this_type& operator=(const this_type& other);

      void set_indices(int l);

      template<typename parameter_type>
      void initialize(parameter_type& parameters);

      void sum_to(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                  FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                  FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void sum_to(FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> >& G2_nonlocal_nu_nu_r_r_w_w,
                  FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> >& G2_nonlocal_nu_nu_k_k_w_w);

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

      nu_r_dmn_type    nu_r_dmn;
      nu_nu_r_dmn_type nu_nu_r_dmn;

      matrix_type tmp;

      matrix_type annihilation_bsr_i;
      matrix_type creation_bsr_j;

      matrix_type overlap_0;
      matrix_type overlap_1;

      FUNC_LIB::function<scalar_type , t>      tau;
      FUNC_LIB::function<complex_type, w     > w_im;
      FUNC_LIB::function<complex_type, w_REAL> w_re;

      FUNC_LIB::function<scalar_type , t>      G_tau;
      FUNC_LIB::function<complex_type, w     > G_w_im;
      FUNC_LIB::function<complex_type, w_REAL> G_w_re;

      FUNC_LIB::function<complex_type, dmn_2<w     , dmn_3<nu_dmn, nu_dmn, r_dmn> > > G_w_im__nu_nu_r;
      FUNC_LIB::function<complex_type, dmn_2<w_REAL, dmn_3<nu_dmn, nu_dmn, r_dmn> > > G_w_re__nu_nu_r;
      FUNC_LIB::function< scalar_type, dmn_2<t     , dmn_3<nu_dmn, nu_dmn, r_dmn> > > G_tau__nu_nu_r;

      FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> > G2_nonlocal_nu_nu_r_r_w_w;
      FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> > G2_nonlocal_nu_nu_k_k_w_w;
    };

    template<typename ed_options>
    sp_Greens_function_data<ed_options>::sp_Greens_function_data()
    {}

    template<typename ed_options>
    sp_Greens_function_data<ed_options>::sp_Greens_function_data(const this_type& /*other*/)
    {}

    template<typename ed_options>
    sp_Greens_function_data<ed_options>::~sp_Greens_function_data()
    {}

    template<typename ed_options>
    template<typename parameter_type>
    void sp_Greens_function_data<ed_options>::initialize(parameter_type& parameters)
    {
      //       cout << "\n\n\t" << __FUNCTION__ << "\n\n";

      complex_type I(0,1);

      scalar_type off_set = parameters.get_real_frequencies_off_set();

      for(int t_i = t::dmn_size()/2; t_i < t::dmn_size(); ++t_i){
        tau(t_i) = t::get_elements()[t_i];
        //      cout << t_i << "\t" << tau(t_i) << "\n";
      }

      for(int w_i=0; w_i<w::dmn_size(); ++w_i){
        w_im(w_i) = I*w::get_elements()[w_i];
        //      cout << w_i << "\t" << w_im(w_i) << "\n";
      }

      for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i){
        w_re(w_i) = w_REAL::get_elements()[w_i] + I*off_set;
        //      cout << w_i << "\t" << w_re(w_i) << "\n";
      }

      G_tau__nu_nu_r = 0;  //G_tau__nu_nu_r.print_fingerprint();

      G_w_im__nu_nu_r = 0; //G_w_im__nu_nu_r.print_fingerprint();
      G_w_re__nu_nu_r = 0; //G_w_re__nu_nu_r.print_fingerprint();

      G2_nonlocal_nu_nu_r_r_w_w = 0;
      G2_nonlocal_nu_nu_k_k_w_w = 0;
    }

    template<typename ed_options>
    void sp_Greens_function_data<ed_options>::set_indices(int l)
    {
      int indices[5];

      nu_nu_r_dmn.linind_2_subind(l, indices);

      b_i = indices[0];
      s_i = indices[1];

      b_j = indices[2];
      s_j = indices[3];

      r_i = r_dmn::parameter_type::origin_index();
      r_j = indices[4];

      delta_r = r_j;

      nu_i = bs_dmn(b_i, s_i);
      nu_j = bs_dmn(b_j, s_j);

      bsr_i = bsr_dmn(b_i, s_i, r_i);
      bsr_j = bsr_dmn(b_j, s_j, r_j);
    }

    template<typename ed_options>
    void sp_Greens_function_data<ed_options>::sum_to(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
                                                     FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      for(int ind=0; ind<nu_nu_r_dmn_type::dmn_size(); ind++)
        {
          set_indices(ind);

          std::cout << ind << "\t" << nu_i << ", " << nu_j << ", " << delta_r << std::endl;

          for(int t_i = t::dmn_size()/2; t_i<t::dmn_size(); ++t_i)
            G_r_t(nu_i, nu_j, delta_r, t_i) += G_tau__nu_nu_r(t_i, ind);

          for(int w_i=0; w_i<w::dmn_size(); ++w_i)
            G_r_w_im(nu_i, nu_j, delta_r, w_i) += G_w_im__nu_nu_r(w_i, ind);

          for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i)
            G_r_w_re(nu_i, nu_j, delta_r, w_i) += G_w_re__nu_nu_r(w_i, ind);
        }

    }

     template<typename ed_options>
    void sp_Greens_function_data<ed_options>::sum_to(FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_r_r_w_w,
                                                     FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_k_k_w_w)
     {
       G_nu_nu_r_r_w_w += G2_nonlocal_nu_nu_r_r_w_w;
       G_nu_nu_k_k_w_w += G2_nonlocal_nu_nu_k_k_w_w;
     }


  }

}

#endif
