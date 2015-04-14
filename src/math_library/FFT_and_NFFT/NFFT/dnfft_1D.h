//-*-C++-*-

#ifndef DNFFT_1D_H
#define DNFFT_1D_H

namespace MATH_ALGORITHMS
{
  namespace NFFT
  {
    /*
    enum NFFT_MODE_NAMES {EXACT, LINEAR, CUBIC};

    template<int max_count, int count>
    struct nfft_atomic_convolution
    {
      template<typename scalar_type>
      inline static void execute_linear(scalar_type* f, scalar_type* M, scalar_type* y)
      {
        f[count] += (M[0+2*count]*y[0]+M[1+2*count]*y[1]);

        nfft_atomic_convolution<max_count, count+1>::execute_linear(f, M, y);
      }

      template<typename scalar_type>
      inline static void execute_cubic(scalar_type* f, scalar_type* M, scalar_type* y)
      {
        f[count] += (M[0+4*count]*y[0]+M[1+4*count]*y[1]+M[2+4*count]*y[2]+M[3+4*count]*y[3]);

        nfft_atomic_convolution<max_count, count+1>::execute_cubic(f, M, y);
      }
    };

    template<int max_count>
    struct nfft_atomic_convolution<max_count, max_count>
    {
      template<typename scalar_type>
      inline static void execute_linear(scalar_type* f, scalar_type* y, scalar_type* M)
      {}

      template<typename scalar_type>
      inline static void execute_cubic(scalar_type* f, scalar_type* y, scalar_type* M)
      {}
    };
    */

    /*! \file cached_auxilery_field_values.h
     *
     *  \author Peter Staar
     *
     *  Contains a templated class over the dimension to represent the cached nfft-plan
     *  this class does only 1 FT, at the end of the accumulation if the error is not measured !!!
     */
    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    class dnfft_1D
    {
    public:

      typedef scalartype                             scalar_type;
      typedef dnfft_1D<scalartype, w_dmn_t, p_dmn_t> this_type;

      const static NFFT_MODE_NAMES DEFAULT_NAME = CUBIC;

      const static int DEFAULT_OVER_SAMPLING   = 8;
      const static int DEFAULT_WINDOW_SAMPLING = 32;

      const static NFFT_MODE_NAMES NAME = DEFAULT_NAME;

      const static int OVER_SAMPLING   = DEFAULT_OVER_SAMPLING;
      const static int WINDOW_SAMPLING = DEFAULT_WINDOW_SAMPLING;


      //typedef gaussian_window_function window_function_t;
      typedef kaiser_bessel_function<1> window_function_t;

      typedef dmn_0<nfft_linear_coefficients_domain>         linear_coefficients_dmn_t;
      typedef dmn_0<nfft_cubic_coefficients_domain>          cubic_coefficients_dmn_t;

      typedef dmn_0<nfft_oversampling_domain   <this_type> > oversampling_dmn_t;
      typedef dmn_0<nfft_window_sampling_domain<this_type> > window_sampling_dmn_t;

      typedef dmn_0<nfft_time_domain<PADDED         , this_type> > padded_time_dmn_t;
      typedef dmn_0<nfft_time_domain<LEFT_ORIENTED  , this_type> > left_oriented_time_dmn_t;
      typedef dmn_0<nfft_time_domain<WINDOW_FUNCTION, this_type> > window_function_time_dmn_t;

      typedef dmn_2<oversampling_dmn_t, window_sampling_dmn_t> convolution_time_dmn_t;

      typedef dmn_2<padded_time_dmn_t       , p_dmn_t> padded_time_p_dmn_t;
      typedef dmn_2<left_oriented_time_dmn_t, p_dmn_t> left_oriented_time_p_dmn_t;

    public:

      dnfft_1D();
      ~dnfft_1D();

      void initialize();

      int get_oversampling_factor();
      int get_window_sampling_factor();

      int  get_maximum_frequency();

      void initialize_domains();
      void initialize_functions();

      void accumulate_at(int  coor, scalartype t_val, scalartype f_val);
      void accumulate_at(int* coor, scalartype t_val, scalartype f_val);

      template<typename other_scalartype>
      void finalize(FUNC_LIB::function<std::complex<other_scalartype>, dmn_2<w_dmn_t, p_dmn_t> >& f_w);

    private:

      void convolute_to_f_tau_exact_test                    (int index, scalartype t_val, scalartype f_val);
      void convolute_to_f_tau_fine_linear_interpolation_test(int index, scalartype t_val, scalartype f_val);
      void convolute_to_f_tau_fine_cubic_interpolation_test (int index, scalartype t_val, scalartype f_val);

      void unroll_linear_interpolation_fast(int N, scalartype* f_tmp_ptr, scalartype* matrix_ptr, scalartype* y_ptr);
      void unroll_cubic_interpolation_fast (int N, scalartype* f_tmp_ptr, scalartype* matrix_ptr, scalartype* y_ptr);
      //void unroll_convolution              (int N, scalartype* f_tau_ptr, scalartype* f_tmp_ptr, scalartype  f_val);

      void fold_time_domain_back();

      template<typename other_scalartype>
      void FT_f_tau_to_f_w(FUNC_LIB::function<std::complex<other_scalartype>, dmn_2<w_dmn_t, p_dmn_t> >& f_w);

    private:

      //       NFFT_MODE_NAMES NAME;

      //       int OVER_SAMPLING;
      //       int WINDOW_SAMPLING;

      double SIGMA_WINDOW_SAMPLING;

      std::vector<int>& integer_wave_vectors;

      p_dmn_t p_dmn_t_obj;

      FUNC_LIB::function<scalartype, padded_time_dmn_t>               tau;
      FUNC_LIB::function<scalartype, window_function_time_dmn_t> fine_tau;

      FUNC_LIB::function<             scalartype , padded_time_p_dmn_t>        f_tau;
      FUNC_LIB::function<             scalartype , oversampling_dmn_t>         f_tmp;

      FUNC_LIB::function<             scalartype , left_oriented_time_p_dmn_t> f_tau_left_oriented;
      FUNC_LIB::function<std::complex<scalartype>, left_oriented_time_p_dmn_t> f_omega;

      FUNC_LIB::function<scalartype, dmn_2<oversampling_dmn_t, window_sampling_dmn_t> > convolution_time_values;
      FUNC_LIB::function<scalartype, dmn_2<oversampling_dmn_t, window_sampling_dmn_t> > window_function;

      FUNC_LIB::function<scalartype, dmn_3<linear_coefficients_dmn_t, oversampling_dmn_t, window_sampling_dmn_t> > linear_convolution_matrices;
      FUNC_LIB::function<scalartype, dmn_3<cubic_coefficients_dmn_t , oversampling_dmn_t, window_sampling_dmn_t> > cubic_convolution_matrices;

      FUNC_LIB::function<scalartype, dmn_3<oversampling_dmn_t, linear_coefficients_dmn_t, window_sampling_dmn_t> > linear_convolution_matrices_2;
      FUNC_LIB::function<scalartype, dmn_3<oversampling_dmn_t, cubic_coefficients_dmn_t , window_sampling_dmn_t> > cubic_convolution_matrices_2;

      FUNC_LIB::function<scalartype, w_dmn_t> phi_wn;
    };

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::dnfft_1D():
      //       NAME(CUBIC),

      //       OVER_SAMPLING  (DEFAULT_OVER_SAMPLING),
      //       WINDOW_SAMPLING(DEFAULT_WINDOW_SAMPLING),

      SIGMA_WINDOW_SAMPLING(2),

      integer_wave_vectors(w_dmn_t::parameter_type::get_integer_wave_vectors()),

      tau("tau"),
      fine_tau("fine_tau"),

      f_tau("f_tau"),
      f_tmp("f_tmp"),

      f_tau_left_oriented("f_tau_left_oriented"),
      f_omega("f_omega"),

      convolution_time_values("convolution_time_values"),
      window_function        ("window_function"),

      linear_convolution_matrices("linear_convolution_matrices"),
      cubic_convolution_matrices ("cubic_convolution_matrices"),

      linear_convolution_matrices_2("linear_convolution_matrices_2"),
      cubic_convolution_matrices_2 ("cubic_convolution_matrices_2"),

      phi_wn("phi_wn")
    {
      initialize_domains();

      initialize_functions();
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::~dnfft_1D()
    {}

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    int dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::get_oversampling_factor()
    {
      return OVER_SAMPLING;
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    int dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::get_window_sampling_factor()
    {
      return WINDOW_SAMPLING;
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    int dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::get_maximum_frequency()
    {
      return w_dmn_t::dmn_size()/2;
      //return *std::max_element(integer_wave_vectors.begin(), integer_wave_vectors.end());
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::initialize()
    {
      f_tau = 0.;
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::initialize_domains()
    {
      oversampling_dmn_t   ::parameter_type::initialize(*this);
      window_sampling_dmn_t::parameter_type::initialize(*this);

      nfft_time_domain<LEFT_ORIENTED         , this_type>::initialize(*this);
      nfft_time_domain<PADDED                , this_type>::initialize(*this);
      nfft_time_domain<WINDOW_FUNCTION       , this_type>::initialize(*this);
      nfft_time_domain<FOLDED_WINDOW_FUNCTION, this_type>::initialize(*this);

      {
        tau                    .reset();
        fine_tau               .reset();
        convolution_time_values.reset();

        assert(tau.size()==padded_time_dmn_t::dmn_size());
        for(int l=0; l<tau.size(); l++)
          tau(l) = padded_time_dmn_t::get_elements()[l];

        assert(fine_tau.size()==window_function_time_dmn_t::dmn_size());
        for(int l=0; l<fine_tau.size(); l++)
          fine_tau(l) = window_function_time_dmn_t::get_elements()[l];

        for(int i=0; i<oversampling_dmn_t::dmn_size(); i++)
          for(int j=0; j<window_sampling_dmn_t::dmn_size(); j++)
            convolution_time_values(i,j) = nfft_time_domain<WINDOW_FUNCTION, this_type>::get_elements()[j+i*window_sampling_dmn_t::dmn_size()];
      }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::initialize_functions()
    {
      {
        window_function_t::n = padded_time_dmn_t::dmn_size();
        window_function_t::m = get_oversampling_factor();

        window_function_t::sigma = SIGMA_WINDOW_SAMPLING;
      }

      {
        f_tmp          .reset();
        f_tau          .reset();
        window_function.reset();

        f_tau_left_oriented.reset();
        f_omega            .reset();

        linear_convolution_matrices.reset();
        cubic_convolution_matrices .reset();

        linear_convolution_matrices_2.reset();
        cubic_convolution_matrices_2 .reset();
	
// 	f_tau                     .print_fingerprint();
//      cubic_convolution_matrices.print_fingerprint();

        int         index = 0;
        scalar_type delta = convolution_time_values(0,1)-convolution_time_values(0,0);
        for(int i=0; i<oversampling_dmn_t::dmn_size(); i++)
          {
            for(int j=0; j<window_sampling_dmn_t::dmn_size(); j++)
              {
                assert(abs(convolution_time_values(i,j)-fine_tau(index))<1.e-6);

                scalar_type tau = convolution_time_values(i,j);

                scalar_type f0 = window_function_t::phi_t(tau);
                scalar_type f1 = window_function_t::phi_t(tau+delta);

                scalar_type df0 = window_function_t::d_phi_t(tau);
                scalar_type df1 = window_function_t::d_phi_t(tau+delta);

                scalar_type a =  f0;
                scalar_type b = df0;

                scalar_type c = -( 3.*f0-3.*f1+2.*df0*delta+df1*delta)/std::pow(delta,2);
                scalar_type d = -(-2.*f0+2.*f1-1.*df0*delta-df1*delta)/std::pow(delta,3);

                window_function(i,j) = f0;

                linear_convolution_matrices(0,i,j) =     f0;
                linear_convolution_matrices(1,i,j) = (f1-f0)/delta;

                cubic_convolution_matrices(0,i,j) = a;
                cubic_convolution_matrices(1,i,j) = b;
                cubic_convolution_matrices(2,i,j) = c;
                cubic_convolution_matrices(3,i,j) = d;

                linear_convolution_matrices_2(i,0,j) =     f0       ;//linear_convolution_matrices(0,i,j);
                linear_convolution_matrices_2(i,1,j) = (f1-f0)/delta;linear_convolution_matrices(1,i,j);

                cubic_convolution_matrices_2(i,0,j) = a;//cubic_convolution_matrices(0,i,j);
                cubic_convolution_matrices_2(i,1,j) = b;//cubic_convolution_matrices(1,i,j);
                cubic_convolution_matrices_2(i,2,j) = c;//cubic_convolution_matrices(2,i,j);
                cubic_convolution_matrices_2(i,3,j) = d;//cubic_convolution_matrices(3,i,j);

                index += 1;
              }
          }
      }

      {
        assert(w_dmn_t::dmn_size() == integer_wave_vectors.size());

        for(int l=0; l<w_dmn_t::dmn_size(); l++)
          phi_wn(l) = window_function_t::phi_wn(integer_wave_vectors[l]);
      }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::accumulate_at(int coor, scalartype t_val, scalartype f_val)
    {
      assert(t_val>-0.5-1.e-6 and t_val<0.5+1.e-6);

      switch(NAME)
        {
        case EXACT:
          convolute_to_f_tau_exact_test(coor, t_val, f_val);
          break;

        case LINEAR:
          convolute_to_f_tau_fine_linear_interpolation_test(coor, t_val, f_val);
          break;

        case CUBIC:
          convolute_to_f_tau_fine_cubic_interpolation_test(coor, t_val, f_val);
          break;

        default:
          throw std::logic_error(__FUNCTION__);
        }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::accumulate_at(int* coor, scalartype t_val, scalartype f_val)
    {
      assert(t_val>-0.5-1.e-6 and t_val<0.5+1.e-6);

      int linind=0;
      p_dmn_t_obj.subind_2_linind(coor, linind);

      switch(NAME)
        {
        case EXACT:
          convolute_to_f_tau_exact_test(linind, t_val, f_val);
          break;

        case LINEAR:
          convolute_to_f_tau_fine_linear_interpolation_test(linind, t_val, f_val);
          break;

        case CUBIC:
          convolute_to_f_tau_fine_cubic_interpolation_test(linind, t_val, f_val);
          break;

        default:
          throw std::logic_error(__FUNCTION__);
        }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    template<typename other_scalartype>
    void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::finalize(FUNC_LIB::function<std::complex<other_scalartype>, dmn_2<w_dmn_t, p_dmn_t> >& f_w)
    {
      fold_time_domain_back();

      FT_f_tau_to_f_w(f_w);
    }

    /*************************************************************
     **                                                         **
     **                                                         **
     **                  private functions                      **
     **                                                         **
     **                                                         **
     *************************************************************/

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::convolute_to_f_tau_exact_test(int index, scalartype t_val, scalartype f_val)
    {
      //cout << __FUNCTION__ <<  endl;
      assert(t_val>-0.5-1.e-6 and t_val<0.5+1.e-6);

      const scalartype T_0           = padded_time_dmn_t::parameter_type::first_element();
      const scalartype one_div_Delta = padded_time_dmn_t::parameter_type::get_one_div_Delta();

      int lambda_0 = (t_val-T_0)*one_div_Delta;

      for(int l=-OVER_SAMPLING; l<=OVER_SAMPLING; l++)
        f_tau(lambda_0+l, index) += f_val*window_function_t::phi_t(tau(lambda_0+l)-t_val);
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::convolute_to_f_tau_fine_linear_interpolation_test(int index, scalartype t_val, scalartype f_val)
    {
      //cout << __FUNCTION__ <<  endl;
      assert(t_val>-0.5-1.e-6 and t_val<0.5+1.e-6);

      const scalartype t_0 = window_function_time_dmn_t::parameter_type::first_element();
      const scalartype T_0 = padded_time_dmn_t         ::parameter_type::first_element();

      const scalartype one_div_delta = padded_time_dmn_t::parameter_type::get_one_div_delta();
      const scalartype one_div_Delta = padded_time_dmn_t::parameter_type::get_one_div_Delta();

      int tau_0 = (t_val-T_0)*one_div_Delta;
      int tau_1 = (tau(tau_0)-t_val-t_0)*one_div_delta;

      assert(tau(tau_0)-1.e-10 < t_val && t_val < tau(tau_0+1)+1.e-10);

      scalartype diff_tau = tau(tau_0)-t_val-fine_tau(tau_1);

      assert(diff_tau>-1.e-6 and diff_tau<padded_time_dmn_t::parameter_type::get_delta());

      scalartype y_ptr[2];

      y_ptr[0] = f_val;
      y_ptr[1] = f_val*diff_tau;

      {
        int tau_index       = tau_0-OVER_SAMPLING;
        int delta_tau_index = tau_1-OVER_SAMPLING*WINDOW_SAMPLING;

        int J =  delta_tau_index   %WINDOW_SAMPLING;
        int I = (delta_tau_index-J)/WINDOW_SAMPLING;
        assert(delta_tau_index == I*WINDOW_SAMPLING+J);

        scalartype* f_tau_ptr  = &f_tau(tau_index, index);
        scalartype* matrix_ptr = &linear_convolution_matrices(0,I,J);

        unroll_linear_interpolation_fast(2*OVER_SAMPLING+1, f_tau_ptr, matrix_ptr, y_ptr);
      }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::convolute_to_f_tau_fine_cubic_interpolation_test(int index, scalartype t_val, scalartype f_val)
    {
      //cout << __FUNCTION__ <<  endl;
      assert(t_val>-0.5-1.e-6 and t_val<0.5+1.e-6);

      const scalartype t_0 = window_function_time_dmn_t::parameter_type::first_element();
      const scalartype T_0 = padded_time_dmn_t         ::parameter_type::first_element();

      const scalartype delta = padded_time_dmn_t::parameter_type::get_delta();
      const scalartype Delta = padded_time_dmn_t::parameter_type::get_Delta();

      const scalartype one_div_delta = padded_time_dmn_t::parameter_type::get_one_div_delta();
      const scalartype one_div_Delta = padded_time_dmn_t::parameter_type::get_one_div_Delta();

      int        tau_0     = (t_val-T_0)*one_div_Delta;
      scalartype t0_val_lb = T_0 + tau_0*Delta;

      assert(tau(tau_0)-1.e-6 < t_val && t_val < tau(tau_0+1)+1.e-6);
      assert(abs(tau(tau_0)-t0_val_lb)<1.e-6);

      //int tau_1 = (tau(tau_0)-t_val-t_0)*one_div_delta;
      int        tau_1     = (t0_val_lb-t_val-t_0)*one_div_delta;
      scalartype t1_val_lb = t_0 + tau_1*delta;

      scalartype diff_tau = t0_val_lb-t_val-t1_val_lb;//fine_tau(tau_1);

      assert(diff_tau>-1.e-6 and diff_tau<padded_time_dmn_t::parameter_type::get_delta());

      scalartype y_ptr[4];

      y_ptr[0] =          f_val;
      y_ptr[1] = y_ptr[0]*diff_tau;
      y_ptr[2] = y_ptr[1]*diff_tau;
      y_ptr[3] = y_ptr[2]*diff_tau;

      {
        int tau_index       = tau_0-OVER_SAMPLING;
        int delta_tau_index = tau_1-OVER_SAMPLING*WINDOW_SAMPLING;

        int J =  delta_tau_index   %WINDOW_SAMPLING;
        int I = (delta_tau_index-J)/WINDOW_SAMPLING;
        assert(delta_tau_index == I*WINDOW_SAMPLING+J);

        scalartype* f_tau_ptr  = &f_tau(tau_index, index);
        scalartype* matrix_ptr = &cubic_convolution_matrices(0,I,J);

        unroll_cubic_interpolation_fast(2*OVER_SAMPLING+1, f_tau_ptr, matrix_ptr, y_ptr);
      }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::unroll_linear_interpolation_fast(int N,
                                                                                         scalartype* f_tau_ptr,
                                                                                         scalartype* matrix_ptr,
                                                                                         scalartype* y_ptr)
    {
      switch(N)
        {
        case 0: return;
        case 1:  nfft_atomic_convolution<1 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 2:  nfft_atomic_convolution<2 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 3:  nfft_atomic_convolution<3 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 4:  nfft_atomic_convolution<4 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 5:  nfft_atomic_convolution<5 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 6:  nfft_atomic_convolution<6 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 7:  nfft_atomic_convolution<7 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 8:  nfft_atomic_convolution<8 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 9:  nfft_atomic_convolution<9 ,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 10: nfft_atomic_convolution<10,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 11: nfft_atomic_convolution<11,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 12: nfft_atomic_convolution<12,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 13: nfft_atomic_convolution<13,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 14: nfft_atomic_convolution<14,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 15: nfft_atomic_convolution<15,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 16: nfft_atomic_convolution<16,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 17: nfft_atomic_convolution<17,0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr); return;

        default:
          /*
            for(int l=0; l<2*OVER_SAMPLING+1; l++)
            f_tau_ptr[l] += (matrix_ptr[0+2*l]*y_ptr[0]+matrix_ptr[1+2*l]*y_ptr[1]);
          */

          throw std::logic_error(__FUNCTION__);
        }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::unroll_cubic_interpolation_fast(int N,
                                                                                        scalartype* f_tau_ptr,
                                                                                        scalartype* matrix_ptr,
                                                                                        scalartype* y_ptr)
    {
      switch(N)
        {
        case 0: return;
        case 1:  nfft_atomic_convolution<1 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 2:  nfft_atomic_convolution<2 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 3:  nfft_atomic_convolution<3 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 4:  nfft_atomic_convolution<4 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 5:  nfft_atomic_convolution<5 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 6:  nfft_atomic_convolution<6 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 7:  nfft_atomic_convolution<7 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 8:  nfft_atomic_convolution<8 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 9:  nfft_atomic_convolution<9 ,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 10: nfft_atomic_convolution<10,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 11: nfft_atomic_convolution<11,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 12: nfft_atomic_convolution<12,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 13: nfft_atomic_convolution<13,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 14: nfft_atomic_convolution<14,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 15: nfft_atomic_convolution<15,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 16: nfft_atomic_convolution<16,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;
        case 17: nfft_atomic_convolution<17,0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr); return;

        default:
          /*
            for(int l=0; l<2*OVER_SAMPLING+1; l++)
            f_tau_ptr[l] += (matrix_ptr[0+4*l]*y_ptr[0]+matrix_ptr[1+4*l]*y_ptr[1]+matrix_ptr[2+4*l]*y_ptr[2]+matrix_ptr[3+4*l]*y_ptr[3]);
          */

          throw std::logic_error(__FUNCTION__);
        }
    }

    /*
      template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
      inline void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::unroll_convolution(int N,
      scalartype* f_tau_ptr,
      scalartype* f_tmp_ptr,
      scalartype  f_val)
      {
      int l=0;

      while(true)
      {
      switch(N-l)
      {
      case 0:
      return;

      case 1:
      f_tau_ptr[0] += f_val*f_tmp_ptr[0];
      return;

      case 2:
      f_tau_ptr[0] += f_val*f_tmp_ptr[0];
      f_tau_ptr[1] += f_val*f_tmp_ptr[1];
      return;

      case 3:
      f_tau_ptr[0] += f_val*f_tmp_ptr[0];
      f_tau_ptr[1] += f_val*f_tmp_ptr[1];
      f_tau_ptr[2] += f_val*f_tmp_ptr[2];
      return;

      default:
      f_tau_ptr[0] += f_val*f_tmp_ptr[0];
      f_tau_ptr[1] += f_val*f_tmp_ptr[1];
      f_tau_ptr[2] += f_val*f_tmp_ptr[2];
      f_tau_ptr[3] += f_val*f_tmp_ptr[3];

      f_tmp_ptr += 4;
      f_tau_ptr += 4;

      l += 4;
      }
      }
      }
    */

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::fold_time_domain_back()
    {
      f_tau_left_oriented = 0;

      int N_padded        = nfft_time_domain<PADDED       , this_type>::get_size();
      int N_left_oriented = nfft_time_domain<LEFT_ORIENTED, this_type>::get_size();

      for(int p_ind=0; p_ind<p_dmn_t::dmn_size(); p_ind++)
        {
          for(int t_ind=0; t_ind<N_padded; t_ind++)
            {
              if(t_ind<2*OVER_SAMPLING)
                {
                  f_tau_left_oriented(t_ind-2*OVER_SAMPLING+N_left_oriented, p_ind) += f_tau(t_ind, p_ind);
                }

              if(t_ind>=2*OVER_SAMPLING and t_ind<N_padded-2*OVER_SAMPLING)
                {
                  f_tau_left_oriented(t_ind-2*OVER_SAMPLING, p_ind)                 += f_tau(t_ind, p_ind);
                }

              if(t_ind>=N_padded-2*OVER_SAMPLING)
                {
                  f_tau_left_oriented(t_ind-2*OVER_SAMPLING-N_left_oriented, p_ind) += f_tau(t_ind, p_ind);
                }
            }
        }
    }

    template<typename scalartype, typename w_dmn_t, typename p_dmn_t>
    template<typename other_scalartype>
    void dnfft_1D<scalartype, w_dmn_t, p_dmn_t>::FT_f_tau_to_f_w(FUNC_LIB::function<std::complex<other_scalartype>, dmn_2<w_dmn_t, p_dmn_t> >& f_w)
    {
      int N = nfft_time_domain<LEFT_ORIENTED, this_type>::get_size();

      double*       f_in  = new double[N];
      fftw_complex* f_out = new fftw_complex[N];

      fftw_plan plan = fftw_plan_dft_r2c_1d(N, f_in, f_out, FFTW_ESTIMATE);

      for(int p_ind=0; p_ind<p_dmn_t::dmn_size(); p_ind++)
        {
          for(int t_ind=0; t_ind<N; t_ind++)
            f_in[t_ind] = f_tau_left_oriented(t_ind, p_ind);

          fftw_execute(plan);

          for(int t_ind=0; t_ind<N/2; t_ind++){
            real(f_omega(t_ind, p_ind)) = -f_out[t_ind][0];
            imag(f_omega(t_ind, p_ind)) =  f_out[t_ind][1];
          }

          for(int t_ind=N/2; t_ind<N; t_ind++){
            real(f_omega(t_ind, p_ind)) = -f_out[N-t_ind][0];
            imag(f_omega(t_ind, p_ind)) = -f_out[N-t_ind][1];
          }
        }

      fftw_destroy_plan(plan);

      std::vector<int> w_indices(0);
      for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++){
        for(int t_ind=0; t_ind<N; t_ind++){
          if(integer_wave_vectors[w_ind]==t_ind or integer_wave_vectors[w_ind]+N==t_ind){
            w_indices.push_back(t_ind);
            break;
          }
        }
      }

      for(int p_ind=0; p_ind<p_dmn_t::dmn_size(); p_ind++)
        for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
          f_w(w_ind, p_ind) = f_omega(w_indices[w_ind], p_ind)/phi_wn(w_ind);

      f_w *= 1./N;

      delete [] f_in;
      delete [] f_out;
    }

  }

}

#endif
