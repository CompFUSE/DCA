//-*-C++-*-

#ifndef DCA_QUADRATURE_INTEGRATION_H
#define DCA_QUADRATURE_INTEGRATION_H

namespace DCA
{

  template<typename parameters_type, typename q_dmn_t>
  class quadrature_integration
  {
#include "type_definitions.h"

    //     typedef typename K_dmn::parameter_type k_cluster_type;
    //     const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

    //     typedef dmn_0<coarsegraining_domain<K_dmn, QUADRATURE_K     > > tet_dmn_type;
    //     typedef dmn_0<coarsegraining_domain<K_dmn, QUADRATURE_ORIGIN> > tet_0_dmn_type;

  public:

    //     quadrature_integration(parameters_type& parameters_ref);
    //     ~quadrature_integration();

    template<typename scalar_type>
    static void quadrature_integration_G_q_w_st(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q);

    template<typename scalar_type>
    static void quadrature_integration_G_q_w_mt(int nr_threads,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q);

    template<typename scalar_type>
    static void quadrature_integration_G_q_t_st(scalar_type beta,
						scalar_type f_val,
						scalar_type t_val,
						FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
						FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q);

    template<typename scalar_type>
    static void quadrature_integration_G_q_t_mt(int nr_threads,
						scalar_type beta,
						scalar_type f_val,
						scalar_type t_val,                                                                                         
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
						FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q);

  private:

    template<typename scalar_type>
    static void* quadrature_integration_G_q_w_mt(void* data);

    template<typename scalar_type>
    static void* quadrature_integration_G_q_t_mt(void* data);

    template<typename scalar_type>
    struct quadrature_integration_functions
    {
      quadrature_integration_functions():
        I_q_ptr(NULL),
        H_q_ptr(NULL),
        S_q_ptr(NULL),
        G_q_ptr(NULL)
      {}

      ~quadrature_integration_functions()
      {}

      scalar_type beta;

      scalar_type f_val;
      scalar_type t_val;

      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >* I_q_ptr;
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >* H_q_ptr;
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >* S_q_ptr;
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >* G_q_ptr;
    };

  };

  template<typename parameters_type, typename q_dmn_t>
  template<typename scalar_type>
  void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_st(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q)
  {
    //cout << __FUNCTION__ << "\t" << q_dmn_t::dmn_size() << "\t" << print_time() << "\n";
    
    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>       G_inv("G_inv", nu::dmn_size());
    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<scalar_type> > geinv_obj(G_inv);

    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++){

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_inv(i,j) = I_q(i,j,q_ind)-H_q(i,j,q_ind)-S_q(i,j,q_ind);

//       if(false)
//         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(G_inv);
//       else
//         LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);

      geinv_obj.execute(G_inv);

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_q(i,j,q_ind) = G_inv(i,j);
    }

  }

  template<typename parameters_type, typename q_dmn_t>
  template<typename scalar_type>
  void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(int nr_threads,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q)
  {
    //cout << __FUNCTION__ << "\t" << q_dmn_t::dmn_size() << "\t" << print_time() << "\n";

    G_q = 0.;

    quadrature_integration_functions<scalar_type> quadrature_integration_functions_obj;

    quadrature_integration_functions_obj.I_q_ptr = &I_q;
    quadrature_integration_functions_obj.H_q_ptr = &H_q;
    quadrature_integration_functions_obj.S_q_ptr = &S_q;
    quadrature_integration_functions_obj.G_q_ptr = &G_q;

    COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> parallelization_obj;

    parallelization_obj.execute(nr_threads, quadrature_integration_G_q_w_mt<scalar_type>, (void*) &quadrature_integration_functions_obj);
  }

  template<typename parameters_type, typename q_dmn_t>
  template<typename scalar_type>
  void* quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(void* void_ptr)
  {
    typedef quadrature_integration_functions<scalar_type> quadrature_functions_type;

    COMP_LIB::posix_data*      data_ptr      = static_cast<COMP_LIB::posix_data     *>(void_ptr);
    quadrature_functions_type* functions_ptr = static_cast<quadrature_functions_type*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q = *(functions_ptr->I_q_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q = *(functions_ptr->H_q_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q = *(functions_ptr->S_q_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q = *(functions_ptr->G_q_ptr);

    q_dmn_t             q_dmn;
    std::pair<int, int> q_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, q_dmn);

    {
      LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>       G_inv("G_inv", nu::dmn_size());
      LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<scalar_type> > geinv_obj(G_inv);

      for(int q_ind=q_bounds.first; q_ind<q_bounds.second; q_ind+=1)
        {
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              G_inv(i,j) = I_q(i,j,q_ind)-H_q(i,j,q_ind)-S_q(i,j,q_ind);

//           if(false)
//             LIN_ALG::GEINV<LIN_ALG::CPU>::execute(G_inv);
//           else
//             LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);

	  geinv_obj.execute(G_inv);

          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              G_q(i,j,q_ind) = G_inv(i,j);
        }
    }

    return 0;
  }

  template<typename parameters_type, typename q_dmn_t>
  template<typename scalar_type>
  void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_st(scalar_type beta,
											 scalar_type f_val,
											 scalar_type t_val,											 
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q)
  {
    //cout << __FUNCTION__ << endl;

    G_q = 0.;

    {
      LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> H_m("H_m", nu::dmn_size());

      LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> L("e_l", nu::dmn_size());
      LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> V("V_l", nu::dmn_size());

      LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> G_t("e_l", nu::dmn_size());

      for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++){

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            H_m(i,j) = H_q(i,j,q_ind)-I_q(i,j,q_ind);

        if(false)
          LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', H_m, L, V);
        else
          LIN_ALG::GEEV<LIN_ALG::CPU>::execute_on_Greens_function_matrix('V', 'U', H_m, L, V);

	/*
	{// test GEEV
	  cout << q_ind << endl;
	  
	  for(int j=0; j<nu::dmn_size(); j++){
	    for(int i=0; i<nu::dmn_size(); i++){

	      std::complex<scalar_type> H_val_0 = H_q(i,j,q_ind)-I_q(i,j,q_ind);

// 	      std::complex<scalar_type> H_val_1 = 0;
// 	      for(int l=0; l<nu::dmn_size(); l++)
// 		H_val_1 += L[l]*real(conj(V(l,i))*V(l,j));

	      std::complex<scalar_type> H_val_1 = 0;
	      for(int l=0; l<nu::dmn_size(); l++)
		H_val_1 += V(i,l)*L[l]*conj(V(j,l));

	      if(abs(H_val_0-H_val_1)>1.e-6) 
		{ 
		  cout << abs(H_val_0-H_val_1) << endl;

		  for(int i=0; i<nu::dmn_size(); i++){
		    for(int j=0; j<nu::dmn_size(); j++){
		      cout << H_q(i,j,q_ind)-I_q(i,j,q_ind) << "\t";
		    }
		    cout << endl;
		  }
		  cout << endl;

		  for(int i=0; i<nu::dmn_size(); i++){
		    for(int j=0; j<nu::dmn_size(); j++){
		      cout << V(i,j) << "\t";
		    }
		    cout << endl;
		  }
		  cout << endl;

		  for(int i=0; i<nu::dmn_size(); i++){
		    for(int j=0; j<nu::dmn_size(); j++){

		      std::complex<scalar_type> H_val_1 = 0;
		      for(int l=0; l<nu::dmn_size(); l++)
			H_val_1 += V(i,l)*L[l]*conj(V(j,l));

		      cout << H_val_1 << "\t";
		    }
		    cout << endl;
		  }
		  cout << endl;

		  assert(false);
		}
	    }
	  }
	}
	*/

        for(int i=0; i<nu::dmn_size(); i++)
	  {

          if(L[i]<0)
            G_t[i] = f_val*std::exp(L[i]*(beta-t_val))/(std::exp(L[i]*beta)+1.);
          else
            G_t[i] = f_val*std::exp(-L[i]*t_val)/(std::exp(-L[i]*beta)+1.);

          if(G_t[i]!=G_t[i]){
            cout << "\n\t warning in compute_G_q_t --> G_t[i] : " << G_t[i] << "\n";
            cout << "\n\tL[i] : " << L[i] << "\n";
            cout << "\n\tbeta : " << beta << "\n";
            cout << "\n\ttau  : " << t_val << "\n";
            cout << "\n\tstd::exp(L[i]*beta)  : "        << std::exp(L[i]*beta)         << "\n";
            cout << "\n\tstd::exp(L[i]*(beta-t_val)) : " << std::exp(L[i]*(beta-t_val)) << "\n";

            throw std::logic_error(__FUNCTION__);
          }
        }

//         for(int j=0; j<nu::dmn_size(); j++)
//           for(int i=0; i<nu::dmn_size(); i++)
//             for(int l=0; l<nu::dmn_size(); l++)
//               G_q(i,j,q_ind) += G_t[l]*real(conj(V(l,i))*V(l,j));

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            for(int l=0; l<nu::dmn_size(); l++)
              G_q(i,j,q_ind) += V(i,l)*G_t[l]*conj(V(j,l));//G_t[l]*real(conj(V(l,i))*V(l,j));
      }
    }
  }

  template<typename parameters_type, typename q_dmn_t>
  template<typename scalar_type>
  void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_mt(int nr_threads,
											 scalar_type beta,
                                                                                         scalar_type f_val, 
											 scalar_type t_val,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                                                         FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q)
  {
    G_q = 0.;

    quadrature_integration_functions<scalar_type> quadrature_integration_functions_obj;

    quadrature_integration_functions_obj.beta = beta;

    quadrature_integration_functions_obj.t_val = t_val;
    quadrature_integration_functions_obj.f_val = f_val;

    quadrature_integration_functions_obj.I_q_ptr = &I_q;
    quadrature_integration_functions_obj.H_q_ptr = &H_q;
    quadrature_integration_functions_obj.G_q_ptr = &G_q;

    COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> parallelization_obj;

    parallelization_obj.execute(nr_threads, quadrature_integration_G_q_t_mt<scalar_type>, (void*) &quadrature_integration_functions_obj);
  }

  template<typename parameters_type, typename q_dmn_t>
  template<typename scalar_type>
  void* quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_mt(void* void_ptr)
  {
    typedef quadrature_integration_functions<scalar_type> quadrature_functions_type;

    COMP_LIB::posix_data*      data_ptr      = static_cast<COMP_LIB::posix_data     *>(void_ptr);
    quadrature_functions_type* functions_ptr = static_cast<quadrature_functions_type*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    double beta = functions_ptr->beta;

    double t_val = functions_ptr->t_val;
    double f_val = functions_ptr->f_val;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q = *(functions_ptr->I_q_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q = *(functions_ptr->H_q_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q = *(functions_ptr->G_q_ptr);

    q_dmn_t             q_dmn;
    std::pair<int, int> q_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, q_dmn);

    {
      LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> H_m("H_m", nu::dmn_size());

      LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> L("e_l", nu::dmn_size());
      LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> V("V_l", nu::dmn_size());

      LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> G_t("e_l", nu::dmn_size());

      for(int q_ind=q_bounds.first; q_ind<q_bounds.second; q_ind+=1)
        {
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              H_m(i,j) = H_q(i,j,q_ind)-I_q(i,j,q_ind);

          if(false)
            LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', H_m, L, V);
          else
            LIN_ALG::GEEV<LIN_ALG::CPU>::execute_on_Greens_function_matrix('V', 'U', H_m, L, V);

          for(int i=0; i<nu::dmn_size(); i++){

            if(L[i]<0)
              G_t[i] = f_val*std::exp(L[i]*(beta-t_val))/(std::exp(L[i]*beta)+1.);
            else
              G_t[i] = f_val*std::exp(-L[i]*t_val)/(std::exp(-L[i]*beta)+1.);

            if(G_t[i]!=G_t[i]){
              cout << "\n\t warning in compute_G_q_t --> G_t[i] : " << G_t[i] << "\n";
              cout << "\n\tL[i] : " << L[i] << "\n";
              cout << "\n\tbeta : " << beta << "\n";
              cout << "\n\ttau  : " << t_val << "\n";
              cout << "\n\tstd::exp(L[i]*beta)  : "        << std::exp(L[i]*beta)         << "\n";
              cout << "\n\tstd::exp(L[i]*(beta-t_val)) : " << std::exp(L[i]*(beta-t_val)) << "\n";

              throw std::logic_error(__FUNCTION__);
            }
          }

//           for(int j=0; j<nu::dmn_size(); j++)
//             for(int i=0; i<nu::dmn_size(); i++)
//               for(int l=0; l<nu::dmn_size(); l++)
//                 G_q(i,j,q_ind) += G_t[l]*real(conj(V(l,i))*V(l,j));

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            for(int l=0; l<nu::dmn_size(); l++)
              G_q(i,j,q_ind) += V(i,l)*G_t[l]*conj(V(j,l));
        }
    }

    return 0;
  }

}

#endif
