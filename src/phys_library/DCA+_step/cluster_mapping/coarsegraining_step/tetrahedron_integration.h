//-*-C++-*-

#ifndef DCA_TETRAHEDRON_INTEGRATION_H
#define DCA_TETRAHEDRON_INTEGRATION_H

namespace DCA
{

  template<typename parameters_type, typename K_dmn>
  class tetrahedron_integration
  {
#include "type_definitions.h"

    typedef typename K_dmn::parameter_type k_cluster_type;
    const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

    typedef dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_K     > > tet_dmn_type;
    typedef dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN> > tet_0_dmn_type;

  public:

    tetrahedron_integration(parameters_type& parameters_ref);
    ~tetrahedron_integration();

  protected:

    template<typename scalar_type>
    void tetrahedron_integration_st(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                    function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                    function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

    template<typename scalar_type>
    void tetrahedron_integration_mt(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                    function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                    function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

    template<typename scalar_type>
    void tetrahedron_integration_st_1D(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                       function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                       function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

//     template<typename scalar_type>
//     void tetrahedron_integration_mt_1D(function<             scalar_type , tet_dmn_type>&                 w_tet,
//                                        function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
//                                        function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

    template<typename scalar_type>
    void tetrahedron_integration_st_2D(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                       function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                       function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

//     template<typename scalar_type>
//     void tetrahedron_integration_mt_2D(function<             scalar_type , tet_dmn_type>&                 w_tet,
//                                        function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
//                                        function<std::complex<scalar_type>, dmn_2<nu, nu> >&              G_int);

    template<typename scalar_type>
    void tetrahedron_integration_st_3D(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                       function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                       function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

//     template<typename scalar_type>
//     void tetrahedron_integration_mt_3D(function<             scalar_type , tet_dmn_type>&                 w_tet,
//                                        function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
//                                        function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int);

    //private:

    template<typename scalar_type>
    static void* tetrahedron_integration_mt_1D(void* data);

    template<typename scalar_type>
    static void* tetrahedron_integration_mt_2D(void* data);

    template<typename scalar_type>
    static void* tetrahedron_integration_mt_3D(void* data);

    template<typename scalar_type>
    struct tetrahedron_integration_functions
    {
      tetrahedron_integration_functions(int N):
	w_tet_ptr(NULL),
	G_tet_ptr(NULL),
	G_int_vec(N)
      {}

      ~tetrahedron_integration_functions()
      {}

      function<             scalar_type , tet_dmn_type>*                 w_tet_ptr;
      function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >* G_tet_ptr;

      std::vector<function<std::complex<scalar_type>, dmn_2<nu, nu> > >  G_int_vec;
    };


  private:

    parameters_type& parameters;
  };

  template<typename parameters_type, typename K_dmn>
  tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration(parameters_type& parameters_ref):
    parameters(parameters_ref)
  {}

  template<typename parameters_type, typename K_dmn>
  tetrahedron_integration<parameters_type, K_dmn>::~tetrahedron_integration()
  {}

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_st(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                                                                   function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                                                                   function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int)
  {
    switch(DIMENSION)
      {
      case 1:
        tetrahedron_integration_st_1D(w_tet, G_tet, G_int);
        break;

      case 2:
        tetrahedron_integration_st_2D(w_tet, G_tet, G_int);
        break;

      case 3:
        tetrahedron_integration_st_3D(w_tet, G_tet, G_int);
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_mt(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                                                                   function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                                                                   function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int)
  {
    int nr_threads = parameters.get_nr_coarsegraining_threads();

    tetrahedron_integration_functions<scalar_type> tetrahedron_integration_functions_obj(nr_threads);

    tetrahedron_integration_functions_obj.w_tet_ptr = &w_tet;
    tetrahedron_integration_functions_obj.G_tet_ptr = &G_tet;

    //tetrahedron_integration_functions_obj.G_int_vec.resize(nr_threads, G_int);

    COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> parallelization_obj;

    switch(DIMENSION)
      {
//       case 1:
//         parallelization_obj.execute(nr_threads, tetrahedron_integration_mt_1D<scalar_type>, (void*) &tetrahedron_integration_functions_obj);
//         break;

      case 2:
	parallelization_obj.execute(nr_threads, tetrahedron_integration_mt_2D<scalar_type>, (void*) &tetrahedron_integration_functions_obj);
	break;

      case 3:
        parallelization_obj.execute(nr_threads, tetrahedron_integration_mt_3D<scalar_type>, (void*) &tetrahedron_integration_functions_obj);
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j) = 0;

    for(int l=0; l<nr_threads; l++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_int(i,j) += tetrahedron_integration_functions_obj.G_int_vec[l](i,j);

  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_st_1D(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                                                                      function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                                                                      function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int)
  {
    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j)=0;

    function<std::complex<scalar_type>, dmn_2<nu, nu> > G_tmp("G_tmp");

    tetrahedron_integration_data<scalar_type> data_obj(nu::dmn_size());

    for(int tet_ind=0; tet_ind<tet_dmn_type::dmn_size(); tet_ind += 4)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2)+w_tet(tet_ind+3);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);

        tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, &G_tmp(0,0), data_obj);

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G_int(i,j) += G_tmp(i,j);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void* tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_mt_1D(void* void_ptr)
  {
    return 0;
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_st_2D(function<             scalar_type , tet_dmn_type>&                 w_tet,
                                                                                      function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                                                                      function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int)
  {
    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j)=0;

    function<std::complex<scalar_type>, dmn_2<nu, nu> > G_tmp("G_tmp");

    for(int tet_ind=0; tet_ind<tet_dmn_type::dmn_size(); tet_ind += 3)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);

        tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2, &G_tmp(0,0));

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G_int(i,j) += G_tmp(i,j);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void* tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_mt_2D(void* void_ptr)
  {
    typedef tetrahedron_integration_functions<scalar_type> tetrahedron_functions_type;

    COMP_LIB::posix_data*       data_ptr      = static_cast<COMP_LIB::posix_data      *>(void_ptr);
    tetrahedron_functions_type* functions_ptr = static_cast<tetrahedron_functions_type*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    function<             scalar_type , tet_dmn_type>&                 w_tet = *(functions_ptr->w_tet_ptr);
    function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet = *(functions_ptr->G_tet_ptr);
    function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int =  (functions_ptr->G_int_vec[id]);

    tet_dmn_type        tet_dmn;
    std::pair<int, int> tet_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, tet_dmn);

    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j)=0;

    function<std::complex<scalar_type>, dmn_2<nu, nu> > G_tmp("G_tmp");

    tetrahedron_integration_data<scalar_type> data_obj(nu::dmn_size());

    for(int tet_ind=tet_bounds.first; tet_ind<tet_bounds.second; tet_ind+=3)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2);

        std::complex<scalar_type>* G_0 = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1 = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2 = &G_tet(0,0,tet_ind+2);

        tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2, &G_tmp(0,0), data_obj);

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G_int(i,j) += G_tmp(i,j);
      }

    return 0;
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_st_3D(function<             scalar_type ,               tet_dmn_type  >& w_tet,
                                                                                      function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet,
                                                                                      function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int)
  {
    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j)=0;

    function<std::complex<scalar_type>, dmn_2<nu, nu> > G_tmp("G_tmp");

    for(int tet_ind=0; tet_ind<tet_dmn_type::dmn_size(); tet_ind += 4)
      {
        scalar_type volume = (w_tet (tet_ind+0)
			      +w_tet(tet_ind+1)
			      +w_tet(tet_ind+2)
			      +w_tet(tet_ind+3));

	{
	  std::vector<double> x(3);
	  std::vector<double> y(3);
	  std::vector<double> z(3);

	  for(int d=0; d<3; d++){
	    x[d] = tet_dmn_type::get_elements()[tet_ind+1][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
	    y[d] = tet_dmn_type::get_elements()[tet_ind+2][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
	    z[d] = tet_dmn_type::get_elements()[tet_ind+3][d]-tet_dmn_type::get_elements()[tet_ind+0][d];	    
	  }

	  if(abs(volume-VECTOR_OPERATIONS::VOLUME(x,y,z))>1.e-6)
	    {
	      cout << tet_ind << "\t" << volume << "\t" << VECTOR_OPERATIONS::VOLUME(x,y,z) << "\n";

	      VECTOR_OPERATIONS::PRINT(x); cout << "\n";
	      VECTOR_OPERATIONS::PRINT(y); cout << "\n";
	      VECTOR_OPERATIONS::PRINT(z); cout << "\n";
	      cout << "\n\n";

	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+0]); cout << "\t" << G_tet(0,0,tet_ind+0) << "\n";
	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+1]); cout << "\t" << G_tet(0,0,tet_ind+1) << "\n";
	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+2]); cout << "\t" << G_tet(0,0,tet_ind+2) << "\n";
	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+3]); cout << "\t" << G_tet(0,0,tet_ind+3) << "\n";
	      cout << "\n\n";

	      assert(false);
	    }
	}

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
        std::complex<scalar_type>* G_3   = &G_tet(0,0,tet_ind+3);

        tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2, G_3, &G_tmp(0,0));

//         for(int j=0; j<nu::dmn_size(); j++)
//           for(int i=0; i<nu::dmn_size(); i++)
// 	    G_tmp(i,j) = (G_tet (i,j,tet_ind+0)
// 			  +G_tet(i,j,tet_ind+1)
// 			  +G_tet(i,j,tet_ind+2)
// 			  +G_tet(i,j,tet_ind+3))*volume/4.;

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G_int(i,j) += G_tmp(i,j);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void* tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_mt_3D(void* void_ptr)
  {
    typedef tetrahedron_integration_functions<scalar_type> tetrahedron_functions_type;

    COMP_LIB::posix_data*       data_ptr      = static_cast<COMP_LIB::posix_data      *>(void_ptr);
    tetrahedron_functions_type* functions_ptr = static_cast<tetrahedron_functions_type*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    function<             scalar_type , tet_dmn_type>&                 w_tet = *(functions_ptr->w_tet_ptr);
    function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet = *(functions_ptr->G_tet_ptr);
    function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int =  (functions_ptr->G_int_vec[id]);

    tet_dmn_type        tet_dmn;
    std::pair<int, int> tet_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, tet_dmn);

    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j)=0;

    function<std::complex<scalar_type>, dmn_2<nu, nu> > G_tmp("G_tmp");

    tetrahedron_integration_data<scalar_type> data_obj(nu::dmn_size());

    for(int tet_ind=tet_bounds.first; tet_ind<tet_bounds.second; tet_ind += 4)    
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2)+w_tet(tet_ind+3);

	/*
	{
	  std::vector<double> x(3);
	  std::vector<double> y(3);
	  std::vector<double> z(3);

	  for(int d=0; d<3; d++){
	    x[d] = tet_dmn_type::get_elements()[tet_ind+1][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
	    y[d] = tet_dmn_type::get_elements()[tet_ind+2][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
	    z[d] = tet_dmn_type::get_elements()[tet_ind+3][d]-tet_dmn_type::get_elements()[tet_ind+0][d];	    
	  }

	  if(abs(volume-VECTOR_OPERATIONS::VOLUME(x,y,z)/6)>1.e-6)
	    {
	      cout << tet_ind << "\t" << volume << "\t" << VECTOR_OPERATIONS::VOLUME(x,y,z) << "\n";

	      VECTOR_OPERATIONS::PRINT(x); cout << "\n";
	      VECTOR_OPERATIONS::PRINT(y); cout << "\n";
	      VECTOR_OPERATIONS::PRINT(z); cout << "\n";
	      cout << "\n\n";

	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+0]); cout << "\t" << G_tet(0,0,tet_ind+0) << "\n";
	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+1]); cout << "\t" << G_tet(0,0,tet_ind+1) << "\n";
	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+2]); cout << "\t" << G_tet(0,0,tet_ind+2) << "\n";
	      VECTOR_OPERATIONS::PRINT(tet_dmn_type::get_elements()[tet_ind+3]); cout << "\t" << G_tet(0,0,tet_ind+3) << "\n";
	      cout << "\n\n";

	      assert(false);
	    }
	}
	*/
       
        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
        std::complex<scalar_type>* G_3   = &G_tet(0,0,tet_ind+3);

        tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2, G_3, &G_tmp(0,0), data_obj);

//         for(int j=0; j<nu::dmn_size(); j++)
//           for(int i=0; i<nu::dmn_size(); i++)
// 	    G_tmp(i,j) = (G_tet(i,j,tet_ind+0)
// 			  +G_tet(i,j,tet_ind+1)
// 			  +G_tet(i,j,tet_ind+2)
// 			  +G_tet(i,j,tet_ind+3))*volume/4.;

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G_int(i,j) += G_tmp(i,j);
      }

    return 0;
  }

}

#endif
