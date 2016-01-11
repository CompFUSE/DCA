//-*-C++-*-

#ifndef DCA_LATTICE_MAP_TP_H
#define DCA_LATTICE_MAP_TP_H

namespace DCA
{
  /*! \ingroup LATTICE-MAPPING
   *
   *  \author Peter Staar
   *  \brief  This class implements the lattice_map.
   *
   */
  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  class lattice_map_tp
  {
#include "type_definitions.h"

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    struct tmp_cluster_domain
    {
      typedef typename target_k_dmn_t::parameter_type::element_type            element_type;
      typedef typename target_k_dmn_t::parameter_type::dmn_specifications_type dmn_specifications_type;
      
      static int get_size()                            { return target_k_dmn_t::dmn_size();     }
      static std::vector<element_type>& get_elements() { return target_k_dmn_t::get_elements(); }
    };

    typedef b                         b_dmn_t;
    typedef w_VERTEX                  w_dmn_t;
    typedef dmn_0<tmp_cluster_domain> tmp_k_dmn_t;

    typedef dmn_4<b_dmn_t, b_dmn_t, tmp_k_dmn_t   , w_dmn_t> tmp_vector_dmn_t;
    typedef dmn_4<b_dmn_t, b_dmn_t, source_k_dmn_t, w_dmn_t> source_vector_dmn_t;
    typedef dmn_4<b_dmn_t, b_dmn_t, target_k_dmn_t, w_dmn_t> target_vector_dmn_t;

  public:

    lattice_map_tp(parameters_type& parameters_ref);
    ~lattice_map_tp();

    template<typename scalartype>
    void execute(FUNC_LIB::function<std::complex<scalartype>, dmn_2<source_vector_dmn_t, source_vector_dmn_t> >& f_source,
                 FUNC_LIB::function<std::complex<scalartype>, dmn_2<target_vector_dmn_t, target_vector_dmn_t> >& f_target);

  private:

    void initialize();

    template<typename k_dmn_t>
    void plot_function(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w> >& f);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    interpolation_tp<parameters_type, source_k_dmn_t, target_k_dmn_t> interpolation_obj;
    deconvolution_tp<parameters_type, source_k_dmn_t, target_k_dmn_t> deconvolution_obj;
  };

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  lattice_map_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::lattice_map_tp(parameters_type& parameters_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    interpolation_obj(parameters),
    deconvolution_obj(parameters)
  {}

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  lattice_map_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::~lattice_map_tp()
  {}

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  template<typename scalartype>
  void lattice_map_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(FUNC_LIB::function<std::complex<scalartype>, dmn_2<source_vector_dmn_t, source_vector_dmn_t> >& f_source,
										FUNC_LIB::function<std::complex<scalartype>, dmn_2<target_vector_dmn_t, target_vector_dmn_t> >& f_target)
  {
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<tmp_vector_dmn_t, tmp_vector_dmn_t> > f_interp("f_interp");

    {
      if(concurrency.id()==0)
	std::cout << "\n\n start tp-interpolation of Gamma \n\n";
      
      interpolation_obj.execute(f_source, f_target);
    }
    
    {
      if(concurrency.id()==0)
	std::cout << "\n\n start tp-deconvolution of Gamma \n\n";

      for(int i=0; i<f_target.size(); i++)
	f_interp(i) = f_target(i);

      deconvolution_obj.execute(f_interp, f_target);
    }
  }

  template<typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
  template<typename k_dmn_t>
  void lattice_map_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::plot_function(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w> >& f)
  {
    std::vector<double> x(0);
    std::vector<double> y(0);

    std::vector<double> z_re(0);
    std::vector<double> z_im(0);

    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
      {
        x.push_back(k_dmn_t::get_elements()[k_ind][0]);
        y.push_back(k_dmn_t::get_elements()[k_ind][1]);
        z_re.push_back(real(f(0,0,k_ind,w::dmn_size()/2)));
        z_im.push_back(imag(f(0,0,k_ind,w::dmn_size()/2)));
      }

    SHOW::heatmap(x,y,z_re);
    SHOW::heatmap(x,y,z_im);
  }

}

#endif
