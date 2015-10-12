//-*-C++-*-

#ifndef DCA_DECONVOLUTION_SP_H
#define DCA_DECONVOLUTION_SP_H

namespace DCA
{
  /*! \ingroup LATTICE-MAPPING
   *
   *  \author Peter Staar
   *  \brief  This class implements the deconvolution in the lattice-mapping.
   *
   */
  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  class deconvolution_sp : public deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>
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

    deconvolution_sp(parameters_type& parameters_ref);
    ~deconvolution_sp();

    void execute(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, source_k_dmn_t, w> >& f_source,
                 FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& Sigma_interp,
                 FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& Sigma_deconv,
                 FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& f_target);

  private:

    void find_shift(FUNC_LIB::function<std::complex<double>, dmn_3<b, b, w> >& shift,
                    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& Sigma_interp);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;
  };

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::deconvolution_sp(parameters_type& parameters_ref):
    deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>(parameters_ref),

    parameters(parameters_ref),
    concurrency(parameters.get_concurrency())
  {}

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::~deconvolution_sp()
  {}

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  void deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, source_k_dmn_t, w> >& f_source,
                                                                                  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& f_interp,
                                                                                  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& f_approx,
                                                                                  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& f_target)
  {
    FUNC_LIB::function<std::complex<double>, dmn_3<b, b, w> > shift;

    find_shift(shift, f_interp);

    typedef dmn_0<dmn<2, int> >    z;
    typedef dmn_5<z, b, b, s, w> p_dmn_t;

    INFERENCE::Richardson_Lucy_deconvolution<parameters_type, target_k_dmn_t, p_dmn_t> RL_obj(parameters);

    FUNC_LIB::function<double, dmn_2<target_k_dmn_t, p_dmn_t> > S_source("S_source");
    FUNC_LIB::function<double, dmn_2<target_k_dmn_t, p_dmn_t> > S_approx("S_approx");
    FUNC_LIB::function<double, dmn_2<target_k_dmn_t, p_dmn_t> > S_target("S_target");

    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
      for(int k_ind=0; k_ind<target_k_dmn_t::dmn_size(); k_ind++){
        for(int j=0; j<b::dmn_size(); j++){
          for(int i=0; i<b::dmn_size(); i++){
            S_source(k_ind, 0, i, j, 0, w_ind) = real(f_interp(i,0,j,0,k_ind,w_ind)) - real(shift(i,j,w_ind));
            S_source(k_ind, 1, i, j, 0, w_ind) = imag(f_interp(i,0,j,0,k_ind,w_ind)) - imag(shift(i,j,w_ind));
            S_source(k_ind, 0, i, j, 1, w_ind) = real(f_interp(i,1,j,1,k_ind,w_ind)) - real(shift(i,j,w_ind));
            S_source(k_ind, 1, i, j, 1, w_ind) = imag(f_interp(i,1,j,1,k_ind,w_ind)) - imag(shift(i,j,w_ind));
          }
        }
      }
    }

    RL_obj.execute(this->T_symmetrized, S_source, S_approx, S_target);

    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
      for(int k_ind=0; k_ind<target_k_dmn_t::dmn_size(); k_ind++){
        for(int j=0; j<b::dmn_size(); j++){
          for(int i=0; i<b::dmn_size(); i++){
            f_approx(i,0,j,0,k_ind,w_ind).real(S_approx(k_ind, 0, i, j, 0, w_ind) + real(shift(i,j,w_ind)));
            f_approx(i,0,j,0,k_ind,w_ind).imag(S_approx(k_ind, 1, i, j, 0, w_ind) + imag(shift(i,j,w_ind)));
            f_approx(i,1,j,1,k_ind,w_ind).real(S_approx(k_ind, 0, i, j, 1, w_ind) + real(shift(i,j,w_ind)));
            f_approx(i,1,j,1,k_ind,w_ind).imag(S_approx(k_ind, 1, i, j, 1, w_ind) + imag(shift(i,j,w_ind)));
          }
        }
      }
    }

    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
      for(int k_ind=0; k_ind<target_k_dmn_t::dmn_size(); k_ind++){
        for(int j=0; j<b::dmn_size(); j++){
          for(int i=0; i<b::dmn_size(); i++){
            f_target(i,0,j,0,k_ind,w_ind).real(S_target(k_ind, 0, i, j, 0, w_ind) + real(shift(i,j,w_ind)));
            f_target(i,0,j,0,k_ind,w_ind).imag(S_target(k_ind, 1, i, j, 0, w_ind) + imag(shift(i,j,w_ind)));
            f_target(i,1,j,1,k_ind,w_ind).real(S_target(k_ind, 0, i, j, 1, w_ind) + real(shift(i,j,w_ind)));
            f_target(i,1,j,1,k_ind,w_ind).imag(S_target(k_ind, 1, i, j, 1, w_ind) + imag(shift(i,j,w_ind)));
          }
        }
      }
    }

  }

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  void deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::find_shift(FUNC_LIB::function<std::complex<double>, dmn_3<b, b, w> >&                      shift,
                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w> >& Sigma_interp)
  {
    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
      for(int j=0; j<b::dmn_size(); j++){
        for(int i=0; i<b::dmn_size(); i++){
          shift(i,j,w_ind).real(0);
          shift(i,j,w_ind).imag(0);
        }
      }
    }

    for(int w_ind=w::dmn_size()/2; w_ind<w::dmn_size(); w_ind++){
      for(int k_ind=0; k_ind<target_k_dmn_t::dmn_size(); k_ind++){
        for(int j=0; j<b::dmn_size(); j++){
          for(int i=0; i<b::dmn_size(); i++){

            if(w_ind<w::dmn_size()/2)
              {
                shift(i,j,w_ind).real(std::min(real(shift(i,j,w_ind)), real(Sigma_interp(i,j,k_ind,w_ind)) ));
                shift(i,j,w_ind).imag(std::min(imag(shift(i,j,w_ind)), imag(Sigma_interp(i,j,k_ind,w_ind)) ));
              }
            else
              {
                shift(i,j,w_ind).real(std::max(real(shift(i,j,w_ind)), real(Sigma_interp(i,j,k_ind,w_ind)) ));
                shift(i,j,w_ind).imag(std::max(imag(shift(i,j,w_ind)), imag(Sigma_interp(i,j,k_ind,w_ind)) ));
              }
          }
        }
      }
    }

    const double factor = 10.;
    
    for(int w_ind=w::dmn_size()/2; w_ind<w::dmn_size(); w_ind++){
      for(int j=0; j<b::dmn_size(); j++){
	for(int i=0; i<b::dmn_size(); i++){
	  shift(i,j,w_ind) *= factor;
	}
      }
    }
    
  }


}

#endif
