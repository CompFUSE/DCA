//-*-C++-*-

#ifndef COMPUTE_BARE_BUBBLE_H
#define COMPUTE_BARE_BUBBLE_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{

  namespace SERIES_EXPANSION
  {

    enum    channel {ph, pp};
    typedef channel channel_type;

    /*!
     * \class  compute_bubble
     *
     * \author Peter Staar
     *
     * \brief  This class implements the computation of the bubble in the particle hole and particle-particle channel.
     */
    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    class compute_bubble
    {

      typedef typename parameters_type::profiler_type    profiler_type;
      typedef typename parameters_type::concurrency_type concurrency_type;

    public:

      typedef FUNC_LIB::function<std::complex<double>, dmn_4<nu , nu, k_dmn_t, w_dmn_t         > > G_function_type;
      typedef FUNC_LIB::function<std::complex<double>, dmn_4<b_b,b_b, k_dmn_t, w_VERTEX_BOSONIC> >   function_type;

    public:

      compute_bubble(parameters_type& parameters_ref);
      ~compute_bubble();

      function_type& get_function();

      void execute_on_lattice(G_function_type& G);
      void execute_on_cluster(G_function_type& G);

      void threaded_execute_on_cluster(G_function_type& G);

      template<IO::FORMAT DATA_FORMAT>
      void write(IO::writer<DATA_FORMAT>& writer);

    private:

      void execute_on_lattice_ph(G_function_type& S);
      void execute_on_cluster_ph(G_function_type& G);

      static void* threaded_execute_on_cluster_ph(void* data);

      void execute_on_lattice_pp(G_function_type& S);
      void execute_on_cluster_pp(G_function_type& G);

      static void* threaded_execute_on_cluster_pp(void* data);

    private:

      parameters_type&  parameters;
      concurrency_type& concurrency;

    protected:

      FUNC_LIB::function<std::complex<double>, dmn_4<b_b,b_b, k_dmn_t, w_VERTEX_BOSONIC> > chi;

    private:

      struct bubble_data
      {
        G_function_type* G_ptr;
        function_type*   chi_ptr;

	concurrency_type* concurrency_ptr;
      };

    };

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::compute_bubble(parameters_type& parameters_ref):
      parameters (parameters_ref),
      concurrency(parameters_ref.get_concurrency())
    {}

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::~compute_bubble()
    {}

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    template<IO::FORMAT DATA_FORMAT>
    void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::write(IO::writer<DATA_FORMAT>& writer)
    {
      writer.execute(chi);
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    typename compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::function_type& compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::get_function()
    {
      return chi;
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::execute_on_cluster(G_function_type& G)
    {
      switch(channel_value)
        {
        case ph:
          execute_on_cluster_ph(G);
          break;

        case pp:
          execute_on_cluster_pp(G);
          break;

        default:
          throw std::logic_error(__FUNCTION__);
        }
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::threaded_execute_on_cluster(G_function_type& G)
    {
      if(concurrency.id()==0)
        std::cout << "\n\n\t\t" << "threaded_execute_on_cluster compute-bubble" << "\n\n";

      profiler_type profiler("threaded_execute_on_cluster compute-bubble", "HTS", __LINE__);

      {
	int nr_threads = parameters.get_nr_HTS_threads();
	
	bubble_data args;
	
	args.G_ptr   = &G;
	args.chi_ptr = &chi;

	args.concurrency_ptr = &concurrency;
	
	COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> pthreads;
	
	switch(channel_value)
	  {
	  case ph:
	    pthreads.execute(nr_threads, threaded_execute_on_cluster_ph, (void*) &args);
	    break;
	    
	  case pp:
	    pthreads.execute(nr_threads, threaded_execute_on_cluster_pp, (void*) &args);
	    break;
	    
	  default:
	    throw std::logic_error(__FUNCTION__);
	  }
      }

      concurrency.sum(chi);

      {
        double factor = -1./(parameters.get_beta()*k_dmn_t::dmn_size());

        chi *= factor;
      }
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::execute_on_cluster_ph(G_function_type& G)
    {
      //cout << __FUNCTION__ << endl;
      if(concurrency.id()==0)
        std::cout << "\n\n\t\t ph-buble \n\n" << std::endl;

      chi.get_name() = "ph-bubble";

      chi = 0.;

      assert(std::fabs(w_VERTEX_BOSONIC::get_elements()[w_VERTEX_BOSONIC::dmn_size()/2])<1.e-6);

      for(int q_ind=0; q_ind<k_dmn_t::dmn_size(); ++q_ind){
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

          int k_plus_q  = k_dmn_t::parameter_type::add(q_ind, k_ind);

          for(int nu_ind=0; nu_ind<w_VERTEX_BOSONIC::dmn_size(); ++nu_ind){

            int nu_c = (nu_ind-w_VERTEX_BOSONIC::dmn_size()/2);

            for(int w_ind=std::fabs(nu_c); w_ind<w_dmn_t::dmn_size()-std::fabs(nu_c); ++w_ind){

              int w_plus_nu = w_ind+nu_c;

              for(int j1=0; j1<b::dmn_size(); ++j1)
                for(int j0=0; j0<b::dmn_size(); ++j0)
                  for(int i1=0; i1<b::dmn_size(); ++i1)
                    for(int i0=0; i0<b::dmn_size(); ++i0)
                      chi(i0, i1, j0, j1, q_ind, nu_ind) += G(i0, j1, k_ind, w_ind)*G(i1,j0, k_plus_q, w_plus_nu);
            }
          }
        }
      }

      double factor = -1./(parameters.get_beta()*k_dmn_t::dmn_size());
      chi *= factor;
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    void* compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::threaded_execute_on_cluster_ph(void* void_ptr)
    {
      COMP_LIB::posix_data* data_ptr   = static_cast<COMP_LIB::posix_data*>(void_ptr);
      bubble_data*          bubble_ptr = static_cast<bubble_data*         >(data_ptr->args);

      G_function_type& G   = *(bubble_ptr->G_ptr);
      function_type&   chi = *(bubble_ptr->chi_ptr);

      concurrency_type& concurrency = *(bubble_ptr->concurrency_ptr);

      k_dmn_t k_dmn;
      std::pair<int, int> k_bounds = concurrency.get_bounds(k_dmn);

      int id         = data_ptr->id;
      int nr_threads = data_ptr->nr_threads;

      k_dmn_t q_dmn;
      std::pair<int, int> q_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, q_dmn);

      for(int q_ind=q_bounds.first; q_ind<q_bounds.second; ++q_ind)
        {
	  double percentage = double(q_ind-q_bounds.first)/double(q_bounds.second-q_bounds.first);

          if(concurrency.id()==0 and id==0 and ( int(100*percentage) % 10==0 ) )
            std::cout << "\t" << int(100*percentage) << " % finished\t" << print_time() << "\n";

          for(int k_ind=k_bounds.first; k_ind<k_bounds.second; ++k_ind)
            {
              int k_plus_q  = k_dmn_t::parameter_type::add(q_ind, k_ind);

              for(int nu_ind=0; nu_ind<w_VERTEX_BOSONIC::dmn_size(); ++nu_ind){

                int nu_c = (nu_ind-w_VERTEX_BOSONIC::dmn_size()/2);

                for(int w_ind=std::fabs(nu_c); w_ind<w_dmn_t::dmn_size()-std::fabs(nu_c); ++w_ind){

                  int w_plus_nu = w_ind+nu_c;

                  for(int j1=0; j1<b::dmn_size(); ++j1)
                    for(int j0=0; j0<b::dmn_size(); ++j0)
                      for(int i1=0; i1<b::dmn_size(); ++i1)
                        for(int i0=0; i0<b::dmn_size(); ++i0)
                          chi(i0, i1, j0, j1, q_ind, nu_ind) += G(i0, j1, k_ind, w_ind)*G(i1,j0, k_plus_q, w_plus_nu);
                }
              }
            }
        }

      return 0;
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::execute_on_cluster_pp(G_function_type& G)
    {
      //cout << __FUNCTION__ << endl;
      if(concurrency.id()==0)
        std::cout << "\n\n\t\t pp-buble \n\n" << std::endl;

      chi.get_name() = "pp-bubble";

      chi = 0.;

      assert(std::fabs(w_VERTEX_BOSONIC::get_elements()[w_VERTEX_BOSONIC::dmn_size()/2])<1.e-6);

      for(int q_ind=0; q_ind<k_dmn_t::dmn_size(); ++q_ind){
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

          int q_minus_k  = k_dmn_t::parameter_type::subtract(k_ind, q_ind);

          for(int nu_ind=0; nu_ind<w_VERTEX_BOSONIC::dmn_size(); ++nu_ind){

            int nu_c = (nu_ind-w_VERTEX_BOSONIC::dmn_size()/2);

            for(int w_ind=std::fabs(nu_c); w_ind<w_dmn_t::dmn_size()-std::fabs(nu_c); ++w_ind){

              int nu_minus_w = nu_c+(w::dmn_size()-1-w_ind);

              for(int j1=0; j1<b::dmn_size(); ++j1)
                for(int j0=0; j0<b::dmn_size(); ++j0)
                  for(int i1=0; i1<b::dmn_size(); ++i1)
                    for(int i0=0; i0<b::dmn_size(); ++i0)
                      chi(i0, i1, j0, j1, q_ind, nu_ind) += G(i0, j0, k_ind, w_ind)*G(i1, j1, q_minus_k, nu_minus_w);
            }
          }
        }
      }

      double factor = -1./(parameters.get_beta()*k_dmn_t::dmn_size());
      chi *= factor;
    }

    template<channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
    void* compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::threaded_execute_on_cluster_pp(void* void_ptr)
    {
      COMP_LIB::posix_data* data_ptr   = static_cast<COMP_LIB::posix_data*>(void_ptr);
      bubble_data*          bubble_ptr = static_cast<bubble_data*         >(data_ptr->args);

      G_function_type& G   = *(bubble_ptr->G_ptr);
      function_type&   chi = *(bubble_ptr->chi_ptr);

      concurrency_type& concurrency = *(bubble_ptr->concurrency_ptr);

      k_dmn_t k_dmn;
      std::pair<int, int> k_bounds = concurrency.get_bounds(k_dmn);

      int id         = data_ptr->id;
      int nr_threads = data_ptr->nr_threads;

      k_dmn_t q_dmn;
      std::pair<int, int> q_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, q_dmn);

      for(int q_ind=q_bounds.first; q_ind<q_bounds.second; ++q_ind)
        {
	  double percentage = double(q_ind-q_bounds.first)/double(q_bounds.second-q_bounds.first);

          if(concurrency.id()==0 and id==0 and ( int(100*percentage) % 10==0 ) )
            std::cout << "\t" << int(100*percentage) << " % finished\t" << print_time() << "\n";

          for(int k_ind=k_bounds.first; k_ind<k_bounds.second; ++k_ind)
            {
              int q_minus_k = k_dmn_t::parameter_type::subtract(k_ind, q_ind);

              for(int nu_ind=0; nu_ind<w_VERTEX_BOSONIC::dmn_size(); ++nu_ind){

                int nu_c = (nu_ind-w_VERTEX_BOSONIC::dmn_size()/2);

                for(int w_ind=std::fabs(nu_c); w_ind<w_dmn_t::dmn_size()-std::fabs(nu_c); ++w_ind){

                  int nu_minus_w = nu_c+(w::dmn_size()-1-w_ind);

                  for(int j1=0; j1<b::dmn_size(); ++j1)
                    for(int j0=0; j0<b::dmn_size(); ++j0)
                      for(int i1=0; i1<b::dmn_size(); ++i1)
                        for(int i0=0; i0<b::dmn_size(); ++i0)
                          chi(i0, i1, j0, j1, q_ind, nu_ind) += G(i0, j0, k_ind, w_ind)*G(i1, j1, q_minus_k, nu_minus_w);
                }
              }
            }
        }

      return 0;
    }


  }

}

#endif
