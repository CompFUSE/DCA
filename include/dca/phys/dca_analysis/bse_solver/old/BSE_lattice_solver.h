//-*-C++-*-

/*
  Modifications compared with BSE_lattice_solver_original.h
  ==========================================================
  BSE_lattice_solver_Peter1:
  1. N_LAMBDAS = 4096; 2. only symmetrize 5 leading eigenvalues
  2a. print out Gamma_lattice

  ==========================================================
  BSE_lattice_solver_Peter2:
  3. comment out characterize_leading_eigenvectors()
  4. symmetrize Gamma_lattice (peoject away odd-frequency part)

  ==========================================================
  BSE_lattice_solver_Peter3:
  5. use sqrt(chi0)*gamma*sqrt(chi0) instead of gamma*chi0

  ==========================================================
  BSE_lattice_solver_Peter4:
  6. use real sqrt(chi0)*gamma*sqrt(chi0) symmetric matrix for particle-particle SC channel
  7. delete some useless routines
  8. always diagonalize full sqrt(chi0)*gamma*sqrt(chi0)
  9. For real symmetric matrix, do not symmetrize_leading_eigenvectors()
  10. print out Gamma_sym
  11. print out sigma_lattice(k,w)

  ==========================================================
  BSE_lattice_solver_Peter5:
  12. only store N_LAMBDA leading eigenvalues and eigenvectors to decrease size of BSE.hdf5.
  Maybe need all eigenpairs for susceptibility
  this way do not need change N_LAMBDA for different cluster sizes for regular DCA, keep using N_LAMBDAS = 10

  ==========================================================
  BSE_lattice_solver_Peter6:
  13. only do lattice mapping for Gamma, chi0 and diagonalization done via python code
  by comment out // BSE_lattice_solver_obj.diagonolize_Gamma_chi_0(Gamma_lattice, chi_0) in BSE_solver.h
  14. move symmetrize Gamma_lattice for even-frequency part for Gamma_sym into compute_Gamma_lattice_3,
  instead of diagonolize_full_Gamma_chi_0
  15. only write Gamma and sigma_lattice data into hdf5
*/

#ifndef BSE_LATTICE_SOLVER_H
#define BSE_LATTICE_SOLVER_H

namespace DCA
{
  /*!
   *  \author Peter Staar
   */
  template<class parameters_type, class MOMS_type>
  class BSE_lattice_solver
  {
#include "type_definitions.h"

    typedef double scalartype;

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    const static int N_LAMBDAS = 10;
    typedef dmn_0<dmn<N_LAMBDAS, int> > lambda_dmn_type;

    const static int N_CUBIC = 3;
    typedef dmn_0<dmn<N_CUBIC, int> > cubic_harmonics_dmn_type;

    typedef dmn_3<b,b,crystal_harmonics_expansion_dmn_t>           chi_vector_dmn_t;

    typedef dmn_4<b,b,k_DCA                             , w_VERTEX> cluster_eigenvector_dmn_t;
    typedef dmn_4<b,b,k_HOST_VERTEX                     , w_VERTEX> lattice_eigenvector_dmn_t;
    typedef dmn_4<b,b,crystal_harmonics_expansion_dmn_t , w_VERTEX> crystal_eigenvector_dmn_t;
    typedef dmn_4<b,b,  cubic_harmonics_dmn_type        , w_VERTEX>   cubic_eigenvector_dmn_t;

    typedef dmn_2<cluster_eigenvector_dmn_t, cluster_eigenvector_dmn_t> DCA_matrix_dmn_t;
    typedef dmn_2<lattice_eigenvector_dmn_t, lattice_eigenvector_dmn_t> HOST_matrix_dmn_t;
    typedef dmn_2<crystal_eigenvector_dmn_t, crystal_eigenvector_dmn_t> crystal_matrix_dmn_t;

  public:

    BSE_lattice_solver(parameters_type& parameters, MOMS_type& MOMS);
    ~BSE_lattice_solver();

    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& get_Gamma_lattice() { return Gamma_lattice; }

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& writer);

    void compute_chi_0_lattice(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0);

    void compute_Gamma_lattice_1(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster);
    //void compute_Gamma_lattice_2(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster);
    void compute_Gamma_lattice_3(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster);

    void diagonolize_Gamma_chi_0(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                 FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

    FUNC_LIB::function<std::complex<scalartype>, lambda_dmn_type>& get_leading_eigenvalues() { return leading_eigenvalues; };
    
  private:

    void initialize();

    void set_chi_0_matrix(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0);

    void diagonolize_full_Gamma_chi_0_real(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                           FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

    void diagonolize_full_Gamma_chi_0(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                      FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

    void record_eigenvalues_and_eigenvectors(LIN_ALG::vector<scalartype, LIN_ALG::CPU>& L,
                                             LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& VR);

    void record_eigenvalues_and_eigenvectors(LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU>& L,
                                             LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VL,
                                             LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VR);

    void diagonolize_folded_Gamma_chi_0(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                        FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

    void record_eigenvalues_and_folded_eigenvectors(LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU>& L,
                                                    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VL,
                                                    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VR);

    void compute_folded_susceptibility(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice,
                                       LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU>& L,
                                       LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VL,
                                       LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VR);

    void symmetrize_leading_eigenvectors();

    void characterize_leading_eigenvectors();

    void print_on_shell();
    void print_on_shell_ppSC();

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    MOMS_type&       MOMS;

    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>                         Gamma_lattice;
    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>                         Gamma_sym;
    FUNC_LIB::function<std::complex<scalartype>, dmn_4<b_b, b_b, k_HOST_VERTEX, w_VERTEX> > chi_0_function;

    FUNC_LIB::function<std::complex<scalartype>, lambda_dmn_type                                   > leading_eigenvalues;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, cubic_eigenvector_dmn_t  > > leading_symmetry_decomposition;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, lattice_eigenvector_dmn_t> > leading_eigenvectors;

    FUNC_LIB::function<scalartype, lambda_dmn_type                                   > leading_eigenvalues_real;
    FUNC_LIB::function<scalartype, dmn_2<lambda_dmn_type, lattice_eigenvector_dmn_t> > leading_eigenvectors_real;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<chi_vector_dmn_t, chi_vector_dmn_t> >         chi_q;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<k_HOST_VERTEX            , crystal_harmonics_expansion_dmn_t> > psi_k;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lattice_eigenvector_dmn_t, crystal_eigenvector_dmn_t        > > crystal_harmonics;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<k_HOST_VERTEX            , cubic_harmonics_dmn_type> > leading_symmetry_functions;
  };

  template<class parameters_type, class MOMS_type>
  BSE_lattice_solver<parameters_type, MOMS_type>::BSE_lattice_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    MOMS(MOMS_ref),

    Gamma_lattice("Gamma_lattice"),
    Gamma_sym("Gamma_sym"),
    chi_0_function("chi-0-function"),

    leading_eigenvalues ("leading-eigenvalues"),
    leading_symmetry_decomposition("leading-symmetry-decomposition"),
    leading_eigenvectors("leading-eigenvectors"),

    leading_eigenvalues_real("leading-eigenvalues-real"),
    leading_eigenvectors_real("leading-eigenvectors-real"),

    chi_q("chi(q)"),

    psi_k("psi_k"),
    crystal_harmonics("crystal_harmonics"),

    leading_symmetry_functions("leading-symmetry-functions")
    //     leading_U_K("leading_U_K"),
    //     leading_Vt_K("leading_Vt_K"),

    //     leading_U_k("leading_U_k_interpolated"),
    //     leading_Vt_k("leading_Vt_k_interpolated")
  {
    initialize();
  }

  template<class parameters_type, class MOMS_type>
  BSE_lattice_solver<parameters_type, MOMS_type>::~BSE_lattice_solver()
  {}

  template<class parameters_type, class MOMS_type>
  template<IO::FORMAT DATA_FORMAT>
  void BSE_lattice_solver<parameters_type, MOMS_type>::write(IO::writer<DATA_FORMAT>& writer)
  {
    if (true)
      {
        writer.execute(leading_eigenvalues);
        writer.execute(leading_eigenvectors);

        writer.execute(leading_symmetry_decomposition);
        writer.execute(leading_symmetry_functions);

        writer.execute(Gamma_lattice);
        writer.execute(chi_0_function);
      }

    else
      {
        writer.execute(leading_eigenvalues_real);
        writer.execute(leading_eigenvectors_real);

        writer.execute(Gamma_sym);
        writer.execute(chi_0_function);
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::initialize()
  {
    {
      double r_cut_off = parameters.get_BSE_cut_off_radius();

      crystal_harmonics_expansion::initialize();

      std::vector<std::vector<double> > r_vecs(0);

      for(int l=0; l<crystal_harmonics_expansion::get_size(); l++)
        if(VECTOR_OPERATIONS::NORM(crystal_harmonics_expansion::get_elements()[l])<r_cut_off)
          r_vecs.push_back(crystal_harmonics_expansion::get_elements()[l]);

      sort(r_vecs.begin(), r_vecs.end(), VECTOR_OPERATIONS::HAS_LARGER_NORM<double>);

      crystal_harmonics_expansion::get_size()     = r_vecs.size();
      crystal_harmonics_expansion::get_elements() = r_vecs;

      if(concurrency.id()==concurrency.last())
        {
          std::cout << "\n\n\t crystal-vectors : \n";
          for(int l=0; l<crystal_harmonics_expansion::get_size(); l++){
            std::cout << "\t" << l << "\t";
            VECTOR_OPERATIONS::PRINT(crystal_harmonics_expansion::get_elements()[l]);
            std::cout << std::endl;
          }
        }
    }

    {
      psi_k            .reset();
      crystal_harmonics.reset();

      const std::complex<double> I(0,1);

      for(int l=0; l<crystal_harmonics_expansion::get_size(); l++)
        {
          std::vector<double> r_vec = crystal_harmonics_expansion::get_elements()[l];

          for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
            psi_k(k_ind, l) = std::exp(I*VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, k_HOST_VERTEX::get_elements()[k_ind]))/std::sqrt(double(k_HOST_VERTEX::dmn_size()));

          for(int w_ind=0; w_ind<w_VERTEX::dmn_size(); w_ind++)
            for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
              for(int m_ind=0; m_ind<b::dmn_size(); m_ind++)
                for(int n_ind=0; n_ind<b::dmn_size(); n_ind++)
                  crystal_harmonics(n_ind, m_ind, k_ind, w_ind, n_ind, m_ind, l, w_ind) = psi_k(k_ind, l);
        }
    }

    {
      leading_symmetry_functions    .reset();
      leading_symmetry_decomposition.reset();

      for(int l=0; l<cubic_harmonics_dmn_type::dmn_size(); l++){

        std::complex<double> norm = 0;

        for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++){

          std::complex<double> value = 0;

          double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
          double ky = k_HOST_VERTEX::get_elements()[k_ind][1];

          switch(l)
            {
            case 0: //s-wave
              value = 1;
              break;

            case 1: //p-wave
              value = cos(kx)+cos(ky);
              break;

            case 2: //d-wave
              value = cos(kx)-cos(ky);
              break;

            default:
              throw std::logic_error(__FUNCTION__);
            }

          norm += conj(value)*value;

          leading_symmetry_functions(k_ind, l) = value;
        }

        for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
          leading_symmetry_functions(k_ind, l) /= std::sqrt(1.e-16+real(norm));
      }
    }
  }

  // Peter's modify May 4
  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::compute_chi_0_lattice(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0)
  {
    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

    typedef DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, LIN_ALG::CPU, parameters_type, MOMS_type> HTS_solver_type;

    typedef DCA::lattice_map_sp   <parameters_type, k_DCA, k_HOST> lattice_map_sp_type;

    typedef DCA::coarsegraining_sp<parameters_type, k_DCA        > coarsegraining_sp_type;
    typedef DCA::coarsegraining_tp<parameters_type, k_HOST_VERTEX> coarsegraining_tp_type;

    lattice_map_sp_type lattice_mapping_obj(parameters);

    MOMS.Sigma_lattice_interpolated  = 0.;
    MOMS.Sigma_lattice_coarsegrained = 0.;

    MOMS.Sigma_lattice = 0.;

    if(parameters.use_interpolated_Self_energy())
      {
        // in case we do the analysis with the DCA+

        if(parameters.use_HTS_approximation())
          {
            coarsegraining_sp_type coarsegraining_sp_obj(parameters);

            MOMS_type MOMS_HTS(parameters);

            MOMS_HTS.H_HOST         = MOMS.H_HOST;
            MOMS_HTS.H_interactions = MOMS.H_interactions;
            HTS_solver_type HTS_solver(parameters, MOMS_HTS);

            lattice_mapping_obj.execute_with_HTS_approximation(MOMS_HTS, HTS_solver, coarsegraining_sp_obj,
                                                               MOMS.Sigma,
                                                               MOMS.Sigma_lattice_interpolated,
                                                               MOMS.Sigma_lattice_coarsegrained,
                                                               MOMS.Sigma_lattice);
          }
        else
          {
            lattice_mapping_obj.execute(MOMS.Sigma, MOMS.Sigma_lattice_interpolated, MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
          }


        {
          coarsegraining_tp_type coarsegraining_tp_obj(parameters);

          coarsegraining_tp_obj.execute(MOMS.H_HOST, MOMS.Sigma_lattice, chi_0_function);

          set_chi_0_matrix(chi_0);
        }
      }
    else
      {
        // in case we do the analysis with the DCA

        {
          coarsegraining_tp_type coarsegraining_tp_obj(parameters);

          coarsegraining_tp_obj.execute(MOMS.H_HOST, MOMS.Sigma, chi_0_function);

          set_chi_0_matrix(chi_0);
        }

        /*
          for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
          for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++){
          for(int i=0; i<nu::dmn_size(); i++){
          for(int j=0; j<nu::dmn_size(); j++){

          int K_ind = -1;//

          MOMS.Sigma_lattice(i, j, k_ind, w_ind) = MOMS.Sigma(i, j, K_ind, w_ind);
          }
          }
          }
          }
        */

      }

    /*
      {
      coarsegraining_tp_type coarsegraining_tp_obj(parameters);

      coarsegraining_tp_obj.execute(MOMS.H_HOST, MOMS.Sigma_lattice, chi_0_function);

      set_chi_0_matrix(chi_0);

      //       scalartype renorm = 1./(parameters.get_beta()*k_HOST_VERTEX::dmn_size());

      //       for(int w_ind=0; w_ind<w_VERTEX::dmn_size(); w_ind++)
      //         for(int K_ind=0; K_ind<k_HOST_VERTEX::dmn_size(); K_ind++)

      //           for(int m2=0; m2<b::dmn_size(); m2++)
      //             for(int n2=0; n2<b::dmn_size(); n2++)

      //               for(int m1=0; m1<b::dmn_size(); m1++)
      //                 for(int n1=0; n1<b::dmn_size(); n1++)
      //                   chi_0(n1,m1,K_ind,w_ind, n2,m2,K_ind,w_ind) = renorm*chi_0_function(n1,m1, n2,m2, K_ind,w_ind);
      }
    */

    // if(concurrency.id()==concurrency.last())
    //   std::cout << "\n\nsymmetrize chi_0_lattice according to the symmetry-group\n" << std::endl;
    // symmetrize::execute(chi_0, parameters.get_q_vector());
  }


  /* original code below:
     template<class parameters_type, class MOMS_type>
     void BSE_lattice_solver<parameters_type, MOMS_type>::compute_chi_0_lattice(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0)
     {
     profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

     if(concurrency.id()==concurrency.last())
     std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

     typedef DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, LIN_ALG::CPU, parameters_type, MOMS_type> HTS_solver_type;

     typedef DCA::lattice_map_sp   <parameters_type, k_DCA, k_HOST> lattice_map_sp_type;

     typedef DCA::coarsegraining_sp<parameters_type, k_DCA        > coarsegraining_sp_type;
     typedef DCA::coarsegraining_tp<parameters_type, k_HOST_VERTEX> coarsegraining_tp_type;

     lattice_map_sp_type lattice_mapping_obj(parameters);

     MOMS.Sigma_lattice_interpolated  = 0.;
     MOMS.Sigma_lattice_coarsegrained = 0.;

     MOMS.Sigma_lattice = 0.;

     if(parameters.use_HTS_approximation())
     {
     coarsegraining_sp_type coarsegraining_sp_obj(parameters);

     MOMS_type MOMS_HTS(parameters);

     MOMS_HTS.H_HOST         = MOMS.H_HOST;
     MOMS_HTS.H_interactions = MOMS.H_interactions;

     HTS_solver_type HTS_solver(parameters, MOMS_HTS);

     lattice_mapping_obj.execute_with_HTS_approximation(MOMS_HTS, HTS_solver, coarsegraining_sp_obj,
     MOMS.Sigma,
     MOMS.Sigma_lattice_interpolated,
     MOMS.Sigma_lattice_coarsegrained,
     MOMS.Sigma_lattice);
     }
     else
     {
     // Need following condition lines for regular DCA; otherwise lattice_mapping_obj.execute do change Sigma
     if (!parameters.use_interpolated_Self_energy()) {
     for (int i=0; i<k_HOST::dmn_size(); i++)
     for (int j=0; j<w::dmn_size(); j++)
     for(int b1=0; b1<b::dmn_size(); b1++)
     for(int b2=0; b2<b::dmn_size(); b2++)
     MOMS.Sigma_lattice(b1,b2,i,j) = MOMS.Sigma(b1,b2,i,j);

     } else {

     lattice_mapping_obj.execute(MOMS.Sigma, MOMS.Sigma_lattice_interpolated, MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
     }
     }

     {
     coarsegraining_tp_type coarsegraining_tp_obj(parameters);

     coarsegraining_tp_obj.execute(MOMS.H_HOST, MOMS.Sigma_lattice, chi_0_function);

     set_chi_0_matrix(chi_0);

     //       scalartype renorm = 1./(parameters.get_beta()*k_HOST_VERTEX::dmn_size());

     //       for(int w_ind=0; w_ind<w_VERTEX::dmn_size(); w_ind++)
     //         for(int K_ind=0; K_ind<k_HOST_VERTEX::dmn_size(); K_ind++)

     //           for(int m2=0; m2<b::dmn_size(); m2++)
     //             for(int n2=0; n2<b::dmn_size(); n2++)

     //               for(int m1=0; m1<b::dmn_size(); m1++)
     //                 for(int n1=0; n1<b::dmn_size(); n1++)
     //                   chi_0(n1,m1,K_ind,w_ind, n2,m2,K_ind,w_ind) = renorm*chi_0_function(n1,m1, n2,m2, K_ind,w_ind);
     }
     }
  */

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::set_chi_0_matrix(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0)
  {
    scalartype renorm = 1./(parameters.get_beta()*k_HOST_VERTEX::dmn_size());

    for(int w_ind=0; w_ind<w_VERTEX::dmn_size(); w_ind++)
      for(int K_ind=0; K_ind<k_HOST_VERTEX::dmn_size(); K_ind++)

        for(int m2=0; m2<b::dmn_size(); m2++)
          for(int n2=0; n2<b::dmn_size(); n2++)

            for(int m1=0; m1<b::dmn_size(); m1++)
              for(int n1=0; n1<b::dmn_size(); n1++)
                chi_0(n1,m1,K_ind,w_ind, n2,m2,K_ind,w_ind) = renorm*chi_0_function(n1,m1, n2,m2, K_ind,w_ind);
  }

  //   template<class parameters_type, class MOMS_type>
  //   void BSE_lattice_solver<parameters_type, MOMS_type>::compute_Gamma_lattice_1(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster)
  //   {
  //     coarsegrain_inversion<parameters_type, k_DCA, k_HOST_VERTEX, QUADRATURE_INTEGRATION> coarsegrain_inversion_obj(parameters);

  //     coarsegrain_inversion_obj.execute(Gamma_cluster, Gamma_lattice, leading_U_K, leading_Vt_K, leading_U_k, leading_Vt_k);
  //   }

  //   template<class parameters_type, class MOMS_type>
  //   void BSE_lattice_solver<parameters_type, MOMS_type>::compute_Gamma_lattice_2(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster)
  //   {
  //     {
  //       if(concurrency.id()==0)
  //         std::cout << "\n\n start tp-interpolation of Gamma \n\n";

  //       interpolation_tp<parameters_type, k_DCA, k_HOST_VERTEX> interpolation_tp_obj(parameters);

  //       interpolation_tp_obj.execute(Gamma_cluster, Gamma_lattice);
  //     }

  //     {
  //       if(concurrency.id()==0)
  //         std::cout << "\n\n start tp-deconvolution of Gamma \n\n";

  //       FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b,b,tmp_cluster_dmn_t,w_VERTEX>, dmn_4<b,b,tmp_cluster_dmn_t,w_VERTEX> > > Gamma_lattice_interp("Gamma_lattice_interp");

  //       for(int i=0; i<Gamma_lattice_interp.size(); i++)
  //         Gamma_lattice_interp(i) = Gamma_lattice(i);

  //       deconvolution_tp<parameters_type, k_DCA, k_HOST_VERTEX> deconvolution_tp_obj(parameters);

  //       deconvolution_tp_obj.execute(Gamma_lattice_interp, Gamma_lattice);
  //     }
  //   }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::compute_Gamma_lattice_3(FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster)
  {
    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    // Need following condition lines for regular DCA; otherwise execute(Gamma_cluster, Gamma_lattice) still do SVD etc.
    if (!parameters.use_interpolated_Self_energy())
      {
        int N = lattice_eigenvector_dmn_t::dmn_size();
        for (int i=0; i<N; i++)
          for (int j=0; j<N; j++)
            Gamma_lattice(i,j) = Gamma_cluster(i,j);
      }
    else
      {
        lattice_map_tp<parameters_type, k_DCA, k_HOST_VERTEX> lattice_map_tp_obj(parameters);

        lattice_map_tp_obj.execute(Gamma_cluster, Gamma_lattice);
      }

    if (false)
    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t symmetrize Gamma_lattice for even-frequency part (only need to symmetrize w2 argument in Gamma(k1,w1,k2,w2)), then compute Gamma_chi_0_lattice " << print_time() << " ...";

      for(int w2=0; w2<w_VERTEX::dmn_size(); w2++)
        for(int K2=0; K2<k_HOST_VERTEX::dmn_size(); K2++)
          for(int m2=0; m2<b::dmn_size(); m2++)
            for(int n2=0; n2<b::dmn_size(); n2++)

              for(int w1=0; w1<w_VERTEX::dmn_size(); w1++)
                for(int K1=0; K1<k_HOST_VERTEX::dmn_size(); K1++)
                  for(int m1=0; m1<b::dmn_size(); m1++)
                    for(int n1=0; n1<b::dmn_size(); n1++)
                    {
                      Gamma_sym(n1,m1,K1,w1,n2,m2,K2,w2) = 0.5* (Gamma_lattice(n1,m1,K1,w1,n2,m2,K2,w2) + Gamma_lattice(n1,m1,K1,w1,n2,m2,K2,w_VERTEX::dmn_size()-1-w2));
                    }
    }
    // if(concurrency.id()==concurrency.last())
    // {
    //     for(int K1=0; K1<k_HOST_VERTEX::dmn_size(); K1++)
    //       for(int w1=0; w1<w_VERTEX::dmn_size(); w1++)
    //       {
    //           std::cout << "\n";
    //           for(int K2=0; K2<k_HOST_VERTEX::dmn_size(); K2++)
    //               for(int w2=0; w2<w_VERTEX::dmn_size(); w2++)
    //               {
    //                   std::cout << real(Gamma_sym(0,0,K1,w1,0,0,K2,w2)) << " ";
    //               }
    //       }
    //     std::cout << "\n";
    // }

    if(parameters.do_symmetrization_of_Gamma())
      {
        if(true)
          {
            concurrency << "symmetrize Gamma_lattice according to the symmetry-group \n\n";

            symmetrize::execute(Gamma_lattice  , parameters.get_q_vector());
        // symmetrize::execute(Gamma_sym, parameters.get_q_vector());
          }

        if(true)
          {
            concurrency << "symmetrize Gamma_lattice according to diagrammatic symmetries \n\n";

            diagrammatic_symmetries<parameters_type> diagrammatic_symmetries_obj(parameters);

            diagrammatic_symmetries_obj.execute(Gamma_lattice);
            // diagrammatic_symmetries_obj.execute(Gamma_symm);
          }
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_Gamma_chi_0(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                                                               FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice)
  {
    if(parameters.do_diagonolization_on_folded_Gamma_chi_0())
      {
        diagonolize_folded_Gamma_chi_0(Gamma_lattice, chi_0_lattice);
      }
    else
      {
        diagonolize_full_Gamma_chi_0(Gamma_lattice, chi_0_lattice);
      }

    characterize_leading_eigenvectors();
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_full_Gamma_chi_0_real(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                                                                         FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice)
  {
    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    int N = lattice_eigenvector_dmn_t::dmn_size();

    LIN_ALG::matrix<scalartype, LIN_ALG::CPU> chi_0_Gamma_chi_0("chi_0_Gamma_chi_0", N);
    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_0_Gamma_chi_0_temp("chi_0_Gamma_chi_0_temp", N);
    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_0_Gamma      ("chi_0_Gamma", N);
    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> Gamma("Gamma", N);
    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_0("chi_0", N);

    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t compute Gamma_chi_0_lattice " << print_time() << " ...";

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          Gamma(i,j) = Gamma_sym(i,j);

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          chi_0(i,j) = sqrt(chi_0_lattice(i,j));

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(chi_0, Gamma, chi_0_Gamma);
      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(chi_0_Gamma, chi_0, chi_0_Gamma_chi_0_temp);

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          chi_0_Gamma_chi_0(i,j) = real(chi_0_Gamma_chi_0_temp(i,j));

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";
    }

    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t diagonalize Gamma_chi_0 in pp SC channel " << print_time() << " ...";

      LIN_ALG::vector<scalartype, LIN_ALG::CPU> L("L (BSE_lattice_solver)", N);

      LIN_ALG::matrix<scalartype, LIN_ALG::CPU> VR("VR (BSE_lattice_solver)", N);

      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', chi_0_Gamma_chi_0, L, VR); // dsyev lapack routine

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";

      record_eigenvalues_and_eigenvectors(L, VR);

      print_on_shell_ppSC();

      //       std::vector<std::pair<std::complex<scalartype>, int> > eigenvals(N);

      //       for(int i=0; i<N; i++){
      //         eigenvals[i].first  = L[i];
      //         eigenvals[i].second = i;
      //       }

      //       stable_sort(eigenvals.begin(), eigenvals.end(), &susceptibility_less_pairs);

      //       { // write down leading eigenvalues ...
      //         for(int i=0; i<N_LAMBDAS; i++)
      //           {
      //             int index = eigenvals[eigenvals.size()-1-i].second;

      //             leading_eigenvalues(i) = L[index];

      //             for(int j=0; j<N; j++)
      //               leading_eigenvectors(i, j) = VR(j, index);
      //           }
      //       }

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";
    }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_full_Gamma_chi_0(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                                                                    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice)
  {
    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    int N = lattice_eigenvector_dmn_t::dmn_size();

    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> Gamma_chi_0("Gamma_chi_0", N);

    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t compute Gamma_chi_0_lattice " << print_time() << " ...";

      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> Gamma("Gamma", N);
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_0("chi_0", N);

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          Gamma(i,j) = Gamma_lattice(i,j);

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          chi_0(i,j) = chi_0_lattice(i,j);

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(Gamma, chi_0, Gamma_chi_0);

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";
    }

    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t diagonalize Gamma_chi_0 " << print_time() << " ...";

      LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU> L("L (BSE_lattice_solver)", N);

      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> VR("VR (BSE_lattice_solver)", N);
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> VL("VL (BSE_lattice_solver)", N);

      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', Gamma_chi_0, L, VL, VR);

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";

      record_eigenvalues_and_eigenvectors(L, VL, VR);

      print_on_shell();

      //       std::vector<std::pair<std::complex<scalartype>, int> > eigenvals(N);

      //       for(int i=0; i<N; i++){
      //         eigenvals[i].first  = L[i];
      //         eigenvals[i].second = i;
      //       }

      //       stable_sort(eigenvals.begin(), eigenvals.end(), &susceptibility_less_pairs);

      //       { // write down leading eigenvalues ...
      //         for(int i=0; i<N_LAMBDAS; i++)
      //           {
      //             int index = eigenvals[eigenvals.size()-1-i].second;

      //             leading_eigenvalues(i) = L[index];

      //             for(int j=0; j<N; j++)
      //               leading_eigenvectors(i, j) = VR(j, index);
      //           }
      //       }

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";
    }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::record_eigenvalues_and_eigenvectors(LIN_ALG::vector<scalartype, LIN_ALG::CPU>& L,
                                                                                           LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& VR)
  {
    int N = lattice_eigenvector_dmn_t::dmn_size();
    //   int M = crystal_eigenvector_dmn_t::dmn_size();

    std::vector<std::pair<scalartype, int> > eigenvals_mod(N);

    for(int i=0; i<N; i++){
      eigenvals_mod[i].first  = std::abs(L[i]-1.);
      eigenvals_mod[i].second = i;
    }

    // sort the eigenvalues by (eig_re-1)**2 + eig_im**2 (ascending order)
    // see src/math_library/static_functions.h for new added real_pair_less
    // replacing original susceptibility_less_pairs
    stable_sort(eigenvals_mod.begin(), eigenvals_mod.end(), &real_pair_less);

    for(int i=0; i<N_LAMBDAS; i++)
      {
        int index = eigenvals_mod[i].second;

        leading_eigenvalues_real(i) = L[index];

        for(int j=0; j<N; j++)
          leading_eigenvectors_real(i, j) = VR(j, index);
      }

    if(concurrency.id()==concurrency.last())
      std::cout << "\n\n\t recording eigenvalues and eigenvectors finished! " << print_time() << "\n";

    //  symmetrize_leading_eigenvectors();
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::record_eigenvalues_and_eigenvectors(LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU>& L,
                                                                                           LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VL,
                                                                                           LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VR)
  {
    int N = lattice_eigenvector_dmn_t::dmn_size();
    int M = crystal_eigenvector_dmn_t::dmn_size();

    std::vector<std::pair<std::complex<scalartype>, int> > eigenvals(M);

    for(int i=0; i<M; i++){
      eigenvals[i].first  = L[i];
      eigenvals[i].second = i;
    }

    stable_sort(eigenvals.begin(), eigenvals.end(), &susceptibility_less_pairs);

    for(int i=0; i<N_LAMBDAS; i++)
    {
      int index = eigenvals[eigenvals.size()-1-i].second;

      leading_eigenvalues(i) = L[index];

      for(int j=0; j<N; j++)
        leading_eigenvectors(i, j) = VR(j, index);
    }

    symmetrize_leading_eigenvectors();
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_folded_Gamma_chi_0(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
                                                                                      FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice)
  {
    if(concurrency.id()==concurrency.last())
      std::cout << __FUNCTION__ << std::endl;

    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    int N = lattice_eigenvector_dmn_t::dmn_size();
    int M = crystal_eigenvector_dmn_t::dmn_size();

    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> Gamma_chi_0_crystal("Gamma_chi_0", M);

    {
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> Gamma_chi_0_lattice("Gamma_chi_0", N);

      {
        if(concurrency.id()==concurrency.last())
          std::cout << "\n\n\t compute Gamma_chi_0_lattice " << print_time() << " ...";

        LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> Gamma("Gamma", N);
        LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_0("chi_0", N);

        for(int j=0; j<N; j++)
          for(int i=0; i<N; i++)
            Gamma(i,j) = Gamma_lattice(i,j);

        for(int j=0; j<N; j++)
          for(int i=0; i<N; i++)
            chi_0(i,j) = chi_0_lattice(i,j);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute(Gamma, chi_0, Gamma_chi_0_lattice);

        if(concurrency.id()==concurrency.last())
          std::cout << " finished " << print_time() << "\n";
      }

      {
        if(concurrency.id()==concurrency.last())
          std::cout << "\n\n\t compute P_Gamma_chi_0_lattice_P " << print_time() << " ...";

        LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> P  ("P"  , std::pair<int,int>(N, M));
        LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> tmp("tmp", std::pair<int,int>(N, M));

        for(int j=0; j<M; j++)
          for(int i=0; i<N; i++)
            P(i,j) = crystal_harmonics(i,j);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'N', Gamma_chi_0_lattice, P, tmp);
        LIN_ALG::GEMM<LIN_ALG::CPU>::execute('C', 'N', P, tmp, Gamma_chi_0_crystal);

        if(concurrency.id()==concurrency.last())
          std::cout << " finished " << print_time() << "\n";
      }
    }

    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t diagonalize P_Gamma_chi_0_lattice_P " << print_time() << " ...";

      LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU> L("L (BSE_lattice_solver)", M);

      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> VL("VL (BSE_lattice_solver)", M);
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> VR("VR (BSE_lattice_solver)", M);

      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', Gamma_chi_0_crystal, L, VL, VR);

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";

      record_eigenvalues_and_folded_eigenvectors(L, VL, VR);

      print_on_shell();

      if(parameters.compute_P_q_lattice()){

      }
    }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::record_eigenvalues_and_folded_eigenvectors(LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU>& L,
                                                                                                  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VL,
                                                                                                  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VR)
  {
    int N = lattice_eigenvector_dmn_t::dmn_size();
    int M = crystal_eigenvector_dmn_t::dmn_size();

    std::vector<std::pair<std::complex<scalartype>, int> > eigenvals(M);

    for(int i=0; i<M; i++){
      eigenvals[i].first  = L[i];
      eigenvals[i].second = i;
    }

    stable_sort(eigenvals.begin(), eigenvals.end(), &susceptibility_less_pairs);

    for(int i=0; i<N_LAMBDAS; i++)
      {
        int index = eigenvals[eigenvals.size()-1-i].second;

        leading_eigenvalues(i) = L[index];

        for(int j=0; j<N; j++)
          for(int l=0; l<M; l++)
            leading_eigenvectors(i, j) += crystal_harmonics(j, l)*VR(l, index);
      }

    symmetrize_leading_eigenvectors();
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::compute_folded_susceptibility(FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice,
                                                                                     LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU>& L,
                                                                                     LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VL,
                                                                                     LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& VR)
  {
    int N = lattice_eigenvector_dmn_t::dmn_size();
    int M = crystal_eigenvector_dmn_t::dmn_size();

    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> P_1_min_Gamma_chi_0_P("P_1_min_Gamma_chi_0_P", M);
    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t invert VR  " << print_time() << " ...";

      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> D_inv   ("D_inv"   , M);
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> VR_D_inv("VR_D_inv", M);

      for(int l=0; l<M; l++)
        D_inv(l,l) = 1./(1.-L[l]);

      VL.copy_from(VR);
      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(VL);

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'N', VR      , D_inv, VR_D_inv);
      LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'N', VR_D_inv, VL   , P_1_min_Gamma_chi_0_P);

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";
    }

    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> P_chi_0_P("P_chi_0_P", M);
    {
      if(concurrency.id()==concurrency.last())
        std::cout << "\n\n\t compute P_chi_0_P " << print_time() << " ...";

      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_0("chi_0", N);

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          chi_0(i,j) = chi_0_lattice(i,j);

      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> P  ("P"  , std::pair<int,int>(N, M));
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> tmp("tmp", std::pair<int,int>(N, M));

      for(int j=0; j<M; j++)
        for(int i=0; i<N; i++)
          P(i,j) = crystal_harmonics(i,j);

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'N', chi_0, P, tmp);
      LIN_ALG::GEMM<LIN_ALG::CPU>::execute('C', 'N', P, tmp, P_chi_0_P);

      if(concurrency.id()==concurrency.last())
        std::cout << " finished " << print_time() << "\n";
    }

    {
      FUNC_LIB::function<std::complex<scalartype>, dmn_2<crystal_eigenvector_dmn_t, crystal_eigenvector_dmn_t> > chi_q_tmp("chi_q_tmp");

      {
        LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> chi_matrix("chi", M);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'N', P_chi_0_P, P_1_min_Gamma_chi_0_P, chi_matrix);

        for(int j=0; j<M; j++)
          for(int i=0; i<M; i++)
            chi_q_tmp(i,j) = chi_matrix(i,j);
      }

      chi_q = 0.;

      for(int w2=0; w2<w_VERTEX::dmn_size(); w2++)
        for(int K2=0; K2<M; K2++)
          for(int m2=0; m2<b::dmn_size(); m2++)
            for(int n2=0; n2<b::dmn_size(); n2++)

              for(int w1=0; w1<w_VERTEX::dmn_size(); w1++)
                for(int K1=0; K1<M; K1++)
                  for(int m1=0; m1<b::dmn_size(); m1++)
                    for(int n1=0; n1<b::dmn_size(); n1++)
                      chi_q(n1,m1,K1, n2,m2,K2) += chi_q_tmp(n1,m1,K1,w1,
                                                             n2,m2,K2,w2);
    }

    if(concurrency.id()==concurrency.last())
      {
        std::cout << "\n\n\t real(chi_q) \n\n";
        for(int i=0; i<M; i++){
          for(int j=0; j<M; j++)
            std::cout << real(chi_q(i,j)) << "\t";
          std::cout << "\n";
        }
        std::cout << "\n";

        std::cout << "\n\n\t imag(chi_q) \n\n";
        for(int i=0; i<M; i++){
          for(int j=0; j<M; j++)
            std::cout << imag(chi_q(i,j)) << "\t";
          std::cout << "\n";
        }
        std::cout << "\n";
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::print_on_shell()
  {
    if(concurrency.id()==concurrency.last())
      std::cout << __FUNCTION__ << std::endl;

    int N = k_HOST_VERTEX::dmn_size();
    int M = crystal_harmonics_expansion_dmn_t::dmn_size();

    if(concurrency.id()==concurrency.last())
      {
        std::cout.precision(6);
        std::cout << std::scientific;

        {
          std::cout << "\n\n\t\t leading eigenvalues : ( T=" << 1./parameters.get_beta() << " )\n\n";
          for(int i=0; i<N_LAMBDAS; i++)
            std::cout << "\t" << i << "\t[" << real(leading_eigenvalues(i)) << ", " << imag(leading_eigenvalues(i)) << "]\n";
        }

        {
          std::cout << "\n\n\t\t leading eigenvectors : \n\n";

          {
            std::cout << "\t" << -1 << "\t[ " << 0. << ", "   << 0. << "]";
            for(int i=0; i<N_LAMBDAS; i++)
              std::cout << "\t[ " << real(leading_eigenvalues(i))
                   << ", "   << imag(leading_eigenvalues(i)) << "]";
            std::cout << "\n\t========================================\n";
          }

          for(int l=0; l<M; l++)
            {
              std::cout << "\t" << l
                   << "\t[ " << crystal_harmonics_expansion::get_elements()[l][0]
                   << ", "   << crystal_harmonics_expansion::get_elements()[l][1] << "]";

              for(int i=0; i<N_LAMBDAS; i++)
                {
                  std::complex<scalartype> result = 0;
                  std::complex<scalartype> norm   = 0;

                  for(int j=0; j<N; j++){
                    result += conj(psi_k(j, l))*leading_eigenvectors(i, 0,0,j,w_VERTEX::dmn_size()/2);
                    norm   += conj(leading_eigenvectors(i, 0,0,j,w_VERTEX::dmn_size()/2))*leading_eigenvectors(i, 0,0,j,w_VERTEX::dmn_size()/2);
                  }

                  std::cout << "\t[ " << real(result/std::sqrt(real(norm)+1.e-16)) << ", " << imag(result/std::sqrt(real(norm)+1.e-16)) << "]";
                }

              std::cout << "\n";
            }
        }
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::print_on_shell_ppSC()
  {
    if(concurrency.id()==concurrency.last())
      std::cout << __FUNCTION__ << std::endl;

    int N = k_HOST_VERTEX::dmn_size();
    //  std::cout << "\t N=" << k_HOST_VERTEX::dmn_size() << "\n";

    if(concurrency.id()==concurrency.last())
      {
        std::cout.precision(6);
        std::cout << std::scientific;
        /*
          {
          std::cout << "\n\n\t\t 10 leading eigenvalues : ( T=" << 1./parameters.get_beta() << " )\n\n";
          for(int i=0; i<10; i++)
          std::cout << "\t" << i << "\t[" << leading_eigenvalues(i) << "]\n";
          }

          {
          std::cout << "\n\n\t\t 10 leading eigenvectors at w=pi*T : \n\n";

          for(int i=0; i<10; i++)
          {
          for(int k_ind=0; k_ind<N; k_ind++)
          {
          double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
          double ky = k_HOST_VERTEX::get_elements()[k_ind][1];
          std::cout << i << "\t" << "k = [" << kx << ",  " << ky << "]\t " << "Phi_k = [" <<
          real(leading_eigenvectors(i,0,0,k_ind,w_VERTEX::dmn_size()/2)) << ",  " <<
          imag(leading_eigenvectors(i,0,0,k_ind,w_VERTEX::dmn_size()/2)) << "]\n";
          }
          std::cout << "----------------------------------------------------------------------" << "\n";
          }
          }
        */
        {
          std::cout << "\n\n\t\t 10 leading eigenvalues : ( T=" << 1./parameters.get_beta() << " )\n\n";
          for(int i=0; i<N_LAMBDAS; i++)
            std::cout << "\t" << i << "\t[" << leading_eigenvalues_real(i) << "]\n";
        }

        {
          int ind0pi = 0;
          int indpi0 = 0;
          // record the index for k=[0,pi] and [pi,0]
          for(int k_ind=0; k_ind<N; k_ind++)
            {
              double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
              double ky = k_HOST_VERTEX::get_elements()[k_ind][1];
              if (std::abs(kx)<1.e-3 && std::abs(ky-3.141592)<1.e-3)
                {
                  ind0pi = k_ind;
                }
              if (std::abs(ky)<1.e-3 && std::abs(kx-3.141592)<1.e-3)
                {
                  indpi0 = k_ind;
                }
            }

          std::cout << "\n\n\t\t Phi_(:,k) for 10 leading eigenvectors: \n\n";
          std::cout << "  w         k=[" << k_HOST_VERTEX::get_elements()[ind0pi][0] << ", " << k_HOST_VERTEX::get_elements()[ind0pi][1] << "]         k=[" << k_HOST_VERTEX::get_elements()[indpi0][0] << ", " << k_HOST_VERTEX::get_elements()[indpi0][1] << "] \n\n";
          for(int i=0; i<N_LAMBDAS; i++)
            {
              for(int w=0; w<w_VERTEX::dmn_size(); w++)
                {
                  std::cout << i << "   " << w_VERTEX::get_elements()[w] << "   " << leading_eigenvectors_real(i,0,0,ind0pi,w) << "   " << leading_eigenvectors_real(i,0,0,indpi0,w) << "\n";
                }
              std::cout << "----------------------------------------------------------------------" << "\n";
            }
        }
        /*
          {
          std::cout << "\n\n\t\t difference between leading eigenvectors at (0,pi) and (pi,0) at w=pi*T : \n\n";
          double eig_0pi=0.;
          double eig_pi0=0.;

          for(int i=0; i<N; i++)
          {
          for(int k_ind=0; k_ind<N; k_ind++)
          {
          double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
          double ky = k_HOST_VERTEX::get_elements()[k_ind][1];

          if ((abs(kx-0.)<0.01)&&(abs(ky-3.14)<0.01)){
          eig_0pi = leading_eigenvectors(w_VERTEX::dmn_size()/2, k_ind, 0,0,i);}
          if ((abs(kx-3.14)<0.01)&&(abs(ky-0.)<0.01)){
          eig_pi0 = leading_eigenvectors(w_VERTEX::dmn_size()/2, k_ind, 0,0,i);}
          }
          std::cout << i << "\t" << eig_0pi-eig_pi0 << "\n";
          }
          }
        */
      }
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::symmetrize_leading_eigenvectors()
  {
    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    if(concurrency.id() == concurrency.first())
      std::cout << "\n\n\t" << __FUNCTION__ << " " << print_time() << std::endl;

    int N = lattice_eigenvector_dmn_t::dmn_size();;

    for(int i=0; i<N_LAMBDAS; i++)
      {
        scalartype N_phases=1.e4;

        std::complex<scalartype> alpha_min=0;
        scalartype norm = 1.e6;

        for(int l=0; l<N_phases; l++)
          {
            std::complex<scalartype> alpha = std::complex<scalartype>(cos(2.*M_PI*l/N_phases), sin(2.*M_PI*l/N_phases));

            scalartype result=0;

            for(int w1=0; w1<w_VERTEX::dmn_size()/2; w1++)
              for(int K1=0; K1<k_HOST_VERTEX::dmn_size(); K1++)
                for(int m1=0; m1<b::dmn_size(); m1++)
                  for(int n1=0; n1<b::dmn_size(); n1++)
                    result += std::abs(alpha*leading_eigenvectors(i, n1,m1,K1,w1)
                                  - conj(alpha*leading_eigenvectors(i, n1,m1,K1,w_VERTEX::dmn_size()-1-w1)));

            if(result < norm){
              norm = result;
              alpha_min = alpha;
            }
          }

        for(int l=0; l<N; l++)
          leading_eigenvectors(i, l) *= alpha_min;
      }

    if(concurrency.id() == concurrency.last())
      std::cout << "\t" << __FUNCTION__ << " finished " << print_time() << std::endl;
  }

  template<class parameters_type, class MOMS_type>
  void BSE_lattice_solver<parameters_type, MOMS_type>::characterize_leading_eigenvectors()
  {
    profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

    for(int n_lam=0; n_lam<N_LAMBDAS; n_lam++)
      {

        for(int w_ind=0; w_ind<w_VERTEX::dmn_size(); w_ind++){
          for(int m_ind=0; m_ind<b::dmn_size(); m_ind++){
            for(int n_ind=0; n_ind<b::dmn_size(); n_ind++){

              for(int n_cub=0; n_cub<N_CUBIC; n_cub++)
                {
                  std::complex<scalartype> scal_prod = 0;

                  std::complex<scalartype> norm_phi  = 0;
                  std::complex<scalartype> norm_psi  = 0;

                  for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++){

                    std::complex<scalartype> phi_k_val = leading_eigenvectors      (n_lam, m_ind, n_ind , k_ind, w_ind);
                    std::complex<scalartype> psi_k_val = leading_symmetry_functions(k_ind, n_cub);

                    scal_prod += conj(psi_k_val)*phi_k_val;

                    norm_phi += conj(phi_k_val)*phi_k_val;
                    norm_psi += conj(psi_k_val)*psi_k_val;
                  }

                  leading_symmetry_decomposition(n_lam, m_ind, n_ind, n_cub, w_ind) = scal_prod/std::sqrt(1.e-16+norm_phi*norm_psi);
                }

            }
          }
        }

      }

  }

}

#endif
