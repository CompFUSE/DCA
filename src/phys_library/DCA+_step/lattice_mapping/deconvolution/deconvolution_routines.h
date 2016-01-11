//-*-C++-*-

#ifndef DCA_DECONVOLUTION_ROUTINES_H
#define DCA_DECONVOLUTION_ROUTINES_H

namespace DCA
{
  /*! \ingroup LATTICE-MAPPING
   *
   *  \author Peter Staar
   *  \brief  This class implements the deconvolution in the lattice-mapping.
   *
   */
  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  class deconvolution_routines
  {
#include "type_definitions.h"

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef typename source_k_dmn_t::parameter_type source_k_cluster_type;
    typedef typename target_k_dmn_t::parameter_type target_k_cluster_type;

    typedef typename source_k_cluster_type::dual_type source_r_cluster_type;
    typedef typename target_k_cluster_type::dual_type target_r_cluster_type;

    typedef dmn_0<source_r_cluster_type>              source_r_dmn_t;
    typedef dmn_0<target_r_cluster_type>              target_r_dmn_t;

    typedef MATH_ALGORITHMS::basis_transform<target_k_dmn_t, target_r_dmn_t> trafo_k_to_r_type;
    typedef MATH_ALGORITHMS::basis_transform<target_r_dmn_t, target_k_dmn_t> trafo_r_to_k_type;

  public:

    deconvolution_routines(parameters_type& parameters_ref);
    ~deconvolution_routines();

    void compute_T_inv_matrix(double epsilon, LIN_ALG::matrix<             double , LIN_ALG::CPU>& T_eps);
    void compute_T_inv_matrix(double epsilon, LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_eps);

  private:

    void initialize();

    void compute_phi_inv(double epsilon);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

  protected:

    FUNC_LIB::function<double, target_r_dmn_t> phi_r;
    FUNC_LIB::function<double, target_r_dmn_t> phi_r_symmetrized;

    FUNC_LIB::function<double, target_r_dmn_t> phi_r_inv;

    LIN_ALG::matrix<double, LIN_ALG::CPU> T;
    LIN_ALG::matrix<double, LIN_ALG::CPU> T_symmetrized;
  };

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::deconvolution_routines(parameters_type& parameters_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    phi_r("phi(r)"),
    phi_r_symmetrized("phi_{sym}(r)"),

    phi_r_inv("phi_r_inv"),

    T            ("T            (deconvolution_routines)", std::pair<int,int>(target_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())),
    T_symmetrized("T_symmetrize (deconvolution_routines)", std::pair<int,int>(target_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size()))
  {
    initialize();
  }

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::~deconvolution_routines()
  {}

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::initialize()
  {
//     if(concurrency.id() == concurrency.first())
//       cout << "\n\n\t initialization of deconvolution_routines has started | time = " << print_time();

    DCA::coarsegraining_sp<parameters_type, source_k_dmn_t> coarsegrain_obj(parameters);

    coarsegrain_obj.compute_phi_r(phi_r);

    phi_r_symmetrized = phi_r;

    symmetrize::execute(phi_r_symmetrized);

    {
      LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_k_to_r = trafo_k_to_r_type::get_transformation_matrix();
      LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_r_to_k = trafo_r_to_k_type::get_transformation_matrix();

      LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> T_k_to_r_scaled("T_k_to_r_scaled");
      LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> T_k_to_k       ("T_k_to_r");

      T_k_to_k.copy_from(T_k_to_r); // resize the matrix;

      {
        T_k_to_r_scaled.copy_from(T_k_to_r);

        for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
          for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
            T_k_to_r_scaled(i,j) *= phi_r(i);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

        for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
          for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
            T(i,j) = real(T_k_to_k(i,j));

        for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
          {
            double result=0;

            for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
              result += T(i,j);

            //    if(abs(result-1)>1.e-6)
            //      cout << "\t T " << i << "\t" << result << "\n";

            //assert(abs(result-1)<1.e-6);
            for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
              T(i,j) /= result;
          }
      }

      //T.print();

      {
        T_k_to_r_scaled.copy_from(T_k_to_r);

        for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
          for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
            T_k_to_r_scaled(i,j) *= phi_r_symmetrized(i);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

        for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
          for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
            T_symmetrized(i,j) = real(T_k_to_k(i,j));

        for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
          {
            double result=0;

            for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
              result += T_symmetrized(i,j);

            //    if(abs(result-1)>1.e-6)
            //      cout << "\t T_symmetrized " << i << "\t" << result << "\n";

            //assert(abs(result-1)<1.e-6);
            for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
              T_symmetrized(i,j) /= result;
          }
      }
    }

    //T_symmetrized.print();

//     if(concurrency.id() == concurrency.first())
//       cout << "\t initialization of deconvolution_routines has ended | time = " << print_time() << "\n";
  }

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_phi_inv(double epsilon)
  {
    if(true)
      {
        for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
          phi_r_inv(i) = std::abs(phi_r(i)) > epsilon ? 1./phi_r(i) : 0.;//1./epsilon;
      }
    else
      {
	for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
	  if(std::abs(phi_r(i)) > epsilon)
	    phi_r_inv(i) = 1./phi_r(i);
	  else
	    phi_r_inv(i) = phi_r(i)>0 ? 1./epsilon : -1./epsilon;
      }
  }

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_T_inv_matrix(double epsilon, LIN_ALG::matrix<double, LIN_ALG::CPU>& T_eps)
  {
    compute_phi_inv(epsilon);

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_k_to_r = trafo_k_to_r_type::get_transformation_matrix();
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_r_to_k = trafo_r_to_k_type::get_transformation_matrix();

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> T_k_to_r_scaled("T_k_to_r_scaled");
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> T_k_to_k       ("T_k_to_r");

    {
      T_k_to_k       .copy_from(T_k_to_r);
      T_k_to_r_scaled.copy_from(T_k_to_r);

      for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
        for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
          T_k_to_r_scaled(i,j) *= phi_r_inv(i);

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

      for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
        for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
          T_eps(i,j) = real(T_k_to_k(i,j));
    }

  }

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_T_inv_matrix(double epsilon, 
												     LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_eps)
  {
    compute_phi_inv(epsilon);

//     FUNC_LIB::function<double, target_r_dmn_t> phi_r_inv;
//     for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
//       phi_r_inv(i) = phi_r(i) > epsilon ? 1./phi_r(i) : 0;//1./epsilon;

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_k_to_r = trafo_k_to_r_type::get_transformation_matrix();
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T_r_to_k = trafo_r_to_k_type::get_transformation_matrix();

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> T_k_to_r_scaled("T_k_to_r_scaled");
//     LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> T_k_to_k       ("T_k_to_r");

    {
//       T_k_to_k       .copy_from(T_k_to_r);
      T_k_to_r_scaled.copy_from(T_k_to_r);

      for(int j=0; j<target_k_dmn_t::dmn_size(); j++)
        for(int i=0; i<target_k_dmn_t::dmn_size(); i++)
          T_k_to_r_scaled(i,j) *= phi_r_inv(i);

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(T_r_to_k, T_k_to_r_scaled, T_eps);
    }

  }

}

#endif
