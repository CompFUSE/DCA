//-*-C++-*-

#ifndef BSE_CLUSTER_SOLVER_H
#define BSE_CLUSTER_SOLVER_H

namespace DCA
{
  /*!
   *  \author Peter Staar
   */
  template<class parameters_type, class MOMS_type>
  class BSE_cluster_solver
  {
#include "type_definitions.h"

    //typedef float scalartype;
    typedef double scalartype;

    typedef typename parameters_type::profiler_type    profiler_t;
    typedef typename parameters_type::concurrency_type concurrency_t;

    typedef dmn_4<b,b,k_DCA,w_VERTEX>                                   cluster_eigenvector_dmn_t;
    typedef dmn_2<cluster_eigenvector_dmn_t, cluster_eigenvector_dmn_t> DCA_matrix_dmn_t;

  public:

    BSE_cluster_solver(parameters_type& parameters, MOMS_type& MOMS);
    ~BSE_cluster_solver();

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);

    void compute_Gamma_cluster();

    FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t             >& get_Gamma_matrix()    { return Gamma_matrix; }
    FUNC_LIB::function<std::complex<double>, dmn_4<b_b, b_b, k_DCA, w_VERTEX> >& get_G_II_0_function() { return G_II_0_function; }

  private:

    void apply_symmetries_sp();
    void apply_symmetries_tp(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
                             FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

    void load_G_II  (FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II);
    void load_G_II_0(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

    void load_G_II_0_function(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

    void solve_BSE_on_cluster(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
			      FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

  private:

    parameters_type& parameters;
    concurrency_t&   concurrency;

    MOMS_type&       MOMS;

    cluster_eigenvector_dmn_t cluster_eigenvector_dmn;

    diagrammatic_symmetries<parameters_type> diagrammatic_symmetries_obj;

    FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t >             Gamma_matrix;
    FUNC_LIB::function<std::complex<double>, dmn_4<b_b, b_b, k_DCA, w_VERTEX> > G_II_0_function;
  };

  template<class parameters_type, class MOMS_type>
  BSE_cluster_solver<parameters_type, MOMS_type>::BSE_cluster_solver(parameters_type& parameters_ref,
								     MOMS_type&       MOMS_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),
    MOMS(MOMS_ref),

    cluster_eigenvector_dmn(),

    diagrammatic_symmetries_obj(parameters),

    Gamma_matrix   ("Gamma_matrix"),
    G_II_0_function("G_II_0_function")
  {}

  template<class parameters_type, class MOMS_type>
  BSE_cluster_solver<parameters_type, MOMS_type>::~BSE_cluster_solver()
  {}

  template<class parameters_type, class MOMS_type>
  template<IO::FORMAT DATA_FORMAT>
  void BSE_cluster_solver<parameters_type, MOMS_type>::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.execute(G_II_0_function);
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::compute_Gamma_cluster()
  {
    FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t> G_II  ("G_II");
    FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t> G_II_0("G_II_0");

    apply_symmetries_sp();

    load_G_II(G_II);

    load_G_II_0         (G_II_0);
    load_G_II_0_function(G_II_0);

    apply_symmetries_tp(G_II, G_II_0);

    solve_BSE_on_cluster(G_II, G_II_0);
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::apply_symmetries_sp()
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << "\n\n";

    profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

    symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

    symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::apply_symmetries_tp(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
                                                                           FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0)
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << "\n\n";

    if(parameters.do_symmetrization_of_Gamma())
      {
        if(true)
          {
            concurrency << "symmetrize Gamma_lattice according to the symmetry-group \n\n";

            symmetrize::execute(G_II  , parameters.get_q_vector());
            symmetrize::execute(G_II_0, parameters.get_q_vector());
          }

        if(true)
          {
            concurrency << "symmetrize Gamma_lattice according to diagrammatic symmetries \n\n";

            diagrammatic_symmetries<parameters_type> diagrammatic_symmetries_obj(parameters);

            diagrammatic_symmetries_obj.execute(G_II  );
            diagrammatic_symmetries_obj.execute(G_II_0);
          }
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::load_G_II(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II)
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << "\n\n";

    int* coor_1 = new int[G_II.signature()];
    int* coor_2 = new int[G_II.signature()];

    for(int i=0; i<MOMS.G4_k_k_w_w.size(); i++)
      {
        MOMS.G4_k_k_w_w.linind_2_subind(i, coor_2);

	// coordinate  0 1 2 3 4 5 6 7
	// G4_k_k_w_w: b b b b k k w w
	// G_II      : b b k w b b k w

        coor_1[0] = coor_2[0];
        coor_1[1] = coor_2[1];
        coor_1[2] = coor_2[4];//k_1
        coor_1[3] = coor_2[6];//w_1
        coor_1[4] = coor_2[2];
        coor_1[5] = coor_2[3];
        coor_1[6] = coor_2[5];//k_2
        coor_1[7] = coor_2[7];//w_2

        G_II(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5], coor_1[6], coor_1[7])
          = MOMS.G4_k_k_w_w(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5], coor_2[6], coor_2[7]);
      }

    delete [] coor_1;
    delete [] coor_2;
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::load_G_II_0(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0)
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << "\n\n";

    G_II_0 = 0.;

    dmn_2<k_DCA, w_VERTEX> k_w_dmn;

    int W        = parameters.get_sp_fermionic_frequencies();
    int W_vertex = w_VERTEX::dmn_size()/2;
    int q        = parameters.get_q_channel();

    int w_nu = parameters.get_w_channel();

    int coor[2];

    for(int i=0; i<k_w_dmn.get_size(); i++)
      {
        k_w_dmn.linind_2_subind(i, coor);

        int k        = coor[0];
        int w_vertex = coor[1];
        int w        = (coor[1] - W_vertex) + W;

        int k_plus_q  = k_DCA::parameter_type::add     (k, q);
        int q_minus_k = k_DCA::parameter_type::subtract(k, q);

        for(int n1=0; n1<b::dmn_size(); n1++){
          for(int n2=0; n2<b::dmn_size(); n2++){
            for(int m1=0; m1<b::dmn_size(); m1++){
              for(int m2=0; m2<b::dmn_size(); m2++){

                switch(parameters.get_vertex_measurement_type())
                  {
                  case PARTICLE_HOLE_TRANSVERSE:
                    {
                      G_II_0(n1,n2,k,w_vertex,
                             m1,m2,k,w_vertex) = -MOMS.G_k_w(n1, e_UP, m2, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w+w_nu);
                      break;
                    }

                  case PARTICLE_HOLE_MAGNETIC:
                    {
                      G_II_0(n1,n2,k,w_vertex,
                             m1,m2,k,w_vertex) = -MOMS.G_k_w(n1, e_UP, m2, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w+w_nu);
                      break;
                    }
                  case PARTICLE_HOLE_CHARGE:
                    {
                      G_II_0(n1,n2,k,w_vertex,
                             m1,m2,k,w_vertex) = -2.*MOMS.G_k_w(n1, e_UP, m1, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m2, e_UP, k_plus_q, w+w_nu);
                      break;
                    }

                  case PARTICLE_PARTICLE_SUPERCONDUCTING:
                    {
                      double wn          = w::get_elements()[w];
                      double w_nu_min_wn = w::get_elements()[w_nu+(2*W-1-w)];
                      double beta        = parameters.get_beta();

                      if(std::fabs( (w_nu*M_PI/beta-wn) - w_nu_min_wn)>1.e-6)
                        throw std::logic_error(__FUNCTION__);


                      G_II_0(n1,n2,k,w_vertex,
                             m1,m2,k,w_vertex) = MOMS.G_k_w(n1, e_UP, m1, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m2, e_UP, q_minus_k, w_nu+(2*W-1-w));
                      break;
                    }

                  default:
                    throw std::logic_error(__FUNCTION__);
                  }
              }
            }
          }
        }
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::load_G_II_0_function(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0)
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << "\n\n";

    for(int w_ind=0; w_ind<w_VERTEX::dmn_size(); w_ind++)
      for(int K_ind=0; K_ind<k_DCA::dmn_size(); K_ind++)

        for(int m2=0; m2<b::dmn_size(); m2++)
          for(int n2=0; n2<b::dmn_size(); n2++)

            for(int m1=0; m1<b::dmn_size(); m1++)
              for(int n1=0; n1<b::dmn_size(); n1++)
                G_II_0_function(n1,m1, n2,m2, K_ind,w_ind) = G_II_0(n1,m1,K_ind,w_ind, n2,m2,K_ind,w_ind);
  }

  template<class parameters_type, class MOMS_type>
  void BSE_cluster_solver<parameters_type, MOMS_type>::solve_BSE_on_cluster(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
                                                                            FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0)
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << "\n\n";

    profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

    scalartype renorm = 1./(parameters.get_beta()*k_DCA::dmn_size());

    int N = cluster_eigenvector_dmn.get_size();

    invert_plan<std::complex<scalartype> > invert_G4  (N);
    invert_plan<std::complex<scalartype> > invert_G4_0(N);

    {
      G_II *= renorm;

      memcpy(&invert_G4.Matrix[0], &G_II(0), sizeof(std::complex<scalartype>)*N*N);
      invert_G4.execute_plan();
    }

    {
      G_II_0 *= renorm;

      memcpy(&invert_G4_0.Matrix[0], &G_II_0(0), sizeof(std::complex<scalartype>)*N*N);
      invert_G4_0.execute_plan();
    }

    for(int j=0; j<N; j++)
      for(int i=0; i<N; i++)
        Gamma_matrix(i, j) = (invert_G4_0.inverted_matrix[i+j*N]
                              -invert_G4 .inverted_matrix[i+j*N]);

    //     if(concurrency.id()==concurrency.last())
    //       std::cout << "symmetrize Gamma_cluster" << std::endl;

    //     symmetrize::execute(Gamma_cluster, MOMS.H_symmetry, parameters.get_q_vector(), false);
    //     diagrammatic_symmetries_obj.execute(Gamma_cluster);
  }

}

#endif
