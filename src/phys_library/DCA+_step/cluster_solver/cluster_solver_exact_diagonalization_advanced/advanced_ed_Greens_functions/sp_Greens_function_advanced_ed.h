//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_SP_GREENS_FUNCTION_H
#define ADVANCED_FERMIONIC_SP_GREENS_FUNCTION_H

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {
    template<typename parameter_type, typename ed_options>
    class fermionic_sp_Greens_function
    {
#include "type_definitions.h"

      //typedef ED_type_definitions<parameter_type, b_dmn, s_dmn, r_dmn> ED_type_def;

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

      typedef typename ed_options::nu_dmn nu_dmn;
      typedef typename ed_options::b_s_r  b_s_r_dmn_type;

      typedef typename ed_options::bs_dmn_type  bs_dmn_type;
      typedef typename ed_options::bsr_dmn_type bsr_dmn_type;

      typedef typename ed_options::nu_nu_r_dmn_type nu_nu_r_dmn_type;

      typedef fermionic_Hamiltonian     <parameter_type, ed_options> fermionic_Hamiltonian_type;
      typedef fermionic_overlap_matrices<parameter_type, ed_options> fermionic_overlap_type;

      typedef Fock_space   <parameter_type, ed_options> fermionic_Fock_space_type;
      typedef Hilbert_space<parameter_type, ed_options> Hilbert_space_type;

      typedef dmn_0<fermionic_Fock_space_type>          fermionic_Fock_dmn_type;

      typedef sp_Greens_function_data<ed_options> sp_Greens_function_data_type;

    public:

      fermionic_sp_Greens_function(parameter_type&             parameters_ref,
                                   fermionic_Hamiltonian_type& Hamiltonian_ref,
                                   fermionic_overlap_type&     overlap_ref);

      ~fermionic_sp_Greens_function();

      template<typename MOMS_w_imag_type, typename MOMS_w_real_type>
      void compute_all_sp_functions(MOMS_w_imag_type& MOMS_imag,
                                    MOMS_w_real_type& MOMS_real,
                                    bool interacting);

      template<typename MOMS_w_imag_type, typename MOMS_w_real_type>
      void compute_all_sp_functions_slow(MOMS_w_imag_type& MOMS_imag,
                                         MOMS_w_real_type& MOMS_real,
                                         bool interacting);

      template<typename w_dmn>
      void compute_S_k_w(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G_k_w,
                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G0_k_w,
                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& S_k_w);

    private:

      /*!
       *  new functions ...
       */
      void compute_real_space_Greens_functions(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
                                               FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
                                               FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                               FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_r_r_w_w,
                                               FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_k_k_w_w);

      void renormalize_real_space_Greens_functions(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
                                                   FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
                                                   FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                                   FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_r_r_w_w,
                                                   FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_k_k_w_w);

      //       template<typename tmp_scalar, typename tmp_dmn>
      //       void compute_momentum_space_Greens_functions(FUNC_LIB::function<tmp_scalar, dmn_4<nu_dmn, nu_dmn, r_dmn, tmp_dmn> >& G_r,
      //                                                    FUNC_LIB::function<tmp_scalar, dmn_4<nu_dmn, nu_dmn, k_dmn, tmp_dmn> >& G_k);

      void compute_Greens_functions_ac_slow(std::vector<sp_Greens_function_data_type>& data_vec);// ,
      //                                             FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
      //                                             FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
      //                                             FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_Greens_functions_ca_slow(std::vector<sp_Greens_function_data_type>& data_vec);// ,
      //                                             FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
      //                                             FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
      //                                             FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_sp_Greens_function(int nu_i_nu_j_delta_r,
                                      scalar_type  E_0,
                                      scalar_type  E_1,
                                      complex_type factor,
                                      sp_Greens_function_data_type& data);


      /*!
       *   Here we try to compute sp-Greens-function as defined by its definition. It might be slower,
       *   but it is a good test to see that we do the right thing for the tp-Greens function
       */
      int has_nonzero_overlap(int HS_i, int HS_j,
                              bool is_creation, int bsr_ind);

      void get_nonzero_overlap(int HS_i, int HS_j,
                               bool is_creation, int bsr_ind,
                               matrix_type& matrix,
                               matrix_type& tmp);

      void compute_sp_permutations(int bsr_0, int bsr_1,
                                   std::vector<std::vector<c_operator> >& sp_perms);

      complex_type phi_e_ikr(int r0, int r1,
                             int k0, int k1);

      complex_type phi_e_iwt(scalar_type E1, scalar_type E2,
                             int         i1, int         i2);

      void compute_sp_Greens_function_slow(int nu_i_nu_j_delta_r,
                                           scalar_type  E_i,
                                           scalar_type  E_j,
                                           complex_type factor,
                                           sp_Greens_function_data_type& data);

      void compute_Greens_functions_slow(std::vector<sp_Greens_function_data_type>& data_vec);

      /*!
       *  old functions ...
       */
      void compute_Greens_functions_main(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                         FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w     > >& G_k_w,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_REAL> >& G_k_w_real,
                                         FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, k_dmn, t     > >& G_k_t);


      void compute_Greens_functions_st   (FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                          FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      //       void compute_Greens_functions_ac_st(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
      //                                           FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
      //                                           FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_Greens_functions_ac_translation(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                   FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                   FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      //       void compute_Greens_functions_ca_st(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
      //                                           FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
      //                                           FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_Greens_functions_ca_translation(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                   FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                   FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_sp_Greens_function(int nu_i, int nu_j, int delta_r,
                                      scalar_type E_0, scalar_type E_1, complex_type factor,
                                      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                      FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);


    private:

      parameter_type&   parameters;
      concurrency_type& concurrency;

      double CUT_OFF;

      fermionic_Hamiltonian_type& Hamiltonian;
      fermionic_overlap_type&     overlap;

      //       std::vector<std::pair<int, int> > ac_Hilbert_space_indices;
      //       std::vector<std::pair<int, int> > ca_Hilbert_space_indices;

      FUNC_LIB::function<vector_type, fermionic_Fock_dmn_type >& eigen_energies;
      FUNC_LIB::function<matrix_type, fermionic_Fock_dmn_type >& eigen_states;

      FUNC_LIB::function<int, dmn_3<fermionic_Fock_dmn_type, fermionic_Fock_dmn_type, b_s_r_dmn_type> >& creation_set_all;
      FUNC_LIB::function<int, dmn_3<fermionic_Fock_dmn_type, fermionic_Fock_dmn_type, b_s_r_dmn_type> >& annihilation_set_all;

      FUNC_LIB::function<int, dmn_2<r_dmn, r_dmn> > rj_minus_ri;

      FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> > G0_nonlocal_nu_nu_r_r_w_w;
      FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> > G0_nonlocal_nu_nu_k_k_w_w;

      FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> > G_nonlocal_nu_nu_r_r_w_w;
      FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> > G_nonlocal_nu_nu_k_k_w_w;
    };

    template<typename parameter_type, typename ed_options>
    fermionic_sp_Greens_function<parameter_type, ed_options>::fermionic_sp_Greens_function(parameter_type&             parameters_ref,
                                                                                           fermionic_Hamiltonian_type& Hamiltonian_ref,
                                                                                           fermionic_overlap_type&     overlap_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CUT_OFF(parameters.get_eigenvalue_cut_off()),

      Hamiltonian(Hamiltonian_ref),
      overlap(overlap_ref),

      //       ac_Hilbert_space_indices(0),
      //       ca_Hilbert_space_indices(0),

      eigen_energies(Hamiltonian.get_eigen_energies()),
      eigen_states  (Hamiltonian.get_eigen_states()),

      creation_set_all    (overlap.get_creation_set_all()),
      annihilation_set_all(overlap.get_annihilation_set_all()),

      rj_minus_ri("rj_minus_ri"),

      G0_nonlocal_nu_nu_r_r_w_w("G0_nonlocal_nu_nu_r_r_w_w"),
      G0_nonlocal_nu_nu_k_k_w_w("G0_nonlocal_nu_nu_k_k_w_w"),

      G_nonlocal_nu_nu_r_r_w_w("G_nonlocal_nu_nu_r_r_w_w"),
      G_nonlocal_nu_nu_k_k_w_w("G_nonlocal_nu_nu_k_k_w_w")
    {
      for(int ri=0; ri<r_dmn::dmn_size(); ri++)
        for(int rj=0; rj<r_dmn::dmn_size(); rj++)
          rj_minus_ri(ri, rj) = r_dmn::parameter_type::subtract(ri, rj);
    }

    template<typename parameter_type, typename ed_options>
    fermionic_sp_Greens_function<parameter_type, ed_options>::~fermionic_sp_Greens_function()
    {}

    //     template<typename parameter_type, typename ed_options>
    //     void fermionic_sp_Greens_function<parameter_type, ed_options>::initialize()
    //     {

    //     }

    template<typename parameter_type, typename ed_options>
    template<typename w_dmn>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_S_k_w(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G_k_w,
                                                                                 FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G0_k_w,
                                                                                 FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& S_k_w)
    {
      if(concurrency.id()==0)
        std::cout << "\n\t" << __FUNCTION__ << std::endl;

      int matrix_size = b::dmn_size()*s::dmn_size()*b::dmn_size()*s::dmn_size();
      int matrix_dim  = b::dmn_size()*s::dmn_size();

      std::complex<double>* G_inverted_matrix                   = new std::complex<double>[matrix_size];
      std::complex<double>* G0_cluster_excluded_inverted_matrix = new std::complex<double>[matrix_size];
      std::complex<double>* Sigma_matrix                        = new std::complex<double>[matrix_size];

      // Sigma = 1/G0 - 1/G

      for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++){
        for(int w_ind=0; w_ind<w_dmn::dmn_size(); w_ind++){

          {
            invert_plan<std::complex<double> > invert_pln(matrix_dim);
            memcpy(invert_pln.Matrix, &G_k_w(0,0,0,0,k_ind,w_ind), sizeof(std::complex<double>)*matrix_size);
            invert_pln.execute_plan();
            memcpy(G_inverted_matrix, invert_pln.inverted_matrix, sizeof(std::complex<double>)*matrix_size);
          }

          {
            invert_plan<std::complex<double> > invert_pln(matrix_dim);
            memcpy(invert_pln.Matrix, &G0_k_w(0,0,0,0,k_ind,w_ind), sizeof(std::complex<double>)*matrix_size);
            invert_pln.execute_plan();
            memcpy(G0_cluster_excluded_inverted_matrix, invert_pln.inverted_matrix, sizeof(std::complex<double>)*matrix_size);
          }

          for(int l=0; l<matrix_size; ++l)
            Sigma_matrix[l] = (G0_cluster_excluded_inverted_matrix[l] - G_inverted_matrix[l]);

          memcpy(&S_k_w(0,0,0,0,k_ind,w_ind), Sigma_matrix, sizeof(std::complex<double>)*matrix_size);
        }
      }

      //SHOW::execute(S_k_w);

      delete [] G_inverted_matrix;
      delete [] G0_cluster_excluded_inverted_matrix;
      delete [] Sigma_matrix;

      if(concurrency.id()==0)
        {
          int N = 8;

//           std::cout << "\n\t G0_k_w\n\n";
//           for(int w_i=w_dmn::dmn_size()/2-N; w_i<w_dmn::dmn_size()/2+N; w_i++){
//             std::cout << w_dmn::get_elements()[w_i] << "\t";
//             for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
//               std::cout << real(G0_k_w(0,0,0,0,k_i,w_i)) << "\t" << imag(G0_k_w(0,0,0,0,k_i,w_i)) << "\t";
//             std::cout << "\n";
//           }
//           std::cout << "\n";

//           std::cout << "\n\t G_k_w\n\n";
//           for(int w_i=w_dmn::dmn_size()/2-N; w_i<w_dmn::dmn_size()/2+N; w_i++){
//             std::cout << w_dmn::get_elements()[w_i] << "\t";
//             for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
//               std::cout << real(G_k_w(0,0,0,0,k_i,w_i)) << "\t" << imag(G_k_w(0,0,0,0,k_i,w_i)) << "\t";
//             std::cout << "\n";
//           }
//           std::cout << "\n";

          std::cout << "\n\t Sigma\n\n";
          for(int w_i=w_dmn::dmn_size()/2-N; w_i<w_dmn::dmn_size()/2+N; w_i++){
            std::cout << w_dmn::get_elements()[w_i] << "\t";
            for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
              std::cout << real(S_k_w(0,0,0,0,k_i,w_i)) << "\t" << imag(S_k_w(0,0,0,0,k_i,w_i)) << "\t";
            std::cout << "\n";
          }
          std::cout << "\n";
        }

      //assert(false);
    }

    template<typename parameter_type, typename ed_options>
    template<typename MOMS_w_imag_type, typename MOMS_w_real_type>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_all_sp_functions_slow(MOMS_w_imag_type& MOMS_imag,
                                                                                                 MOMS_w_real_type& MOMS_real,
                                                                                                 bool interacting)
    {
      if(interacting)
        {
          compute_real_space_Greens_functions    (MOMS_imag.G_r_w, MOMS_real.G_r_w, MOMS_imag.G_r_t, G_nonlocal_nu_nu_r_r_w_w, G_nonlocal_nu_nu_k_k_w_w);
          renormalize_real_space_Greens_functions(MOMS_imag.G_r_w, MOMS_real.G_r_w, MOMS_imag.G_r_t, G_nonlocal_nu_nu_r_r_w_w, G_nonlocal_nu_nu_k_k_w_w);

          //SHOW::execute(MOMS_imag.G_r_t);

          MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(MOMS_imag.G_r_w, MOMS_imag.G_k_w);
          MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(MOMS_real.G_r_w, MOMS_real.G_k_w);
          MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(MOMS_imag.G_r_t, MOMS_imag.G_k_t);

          //SHOW::execute(MOMS_imag.G_r_w);
          //SHOW::execute(MOMS_imag.G_k_w);
        }
      else
        {
          compute_real_space_Greens_functions    (MOMS_imag.G0_r_w, MOMS_real.G0_r_w, MOMS_imag.G0_r_t, G0_nonlocal_nu_nu_r_r_w_w, G0_nonlocal_nu_nu_k_k_w_w);
          renormalize_real_space_Greens_functions(MOMS_imag.G0_r_w, MOMS_real.G0_r_w, MOMS_imag.G0_r_t, G0_nonlocal_nu_nu_r_r_w_w, G0_nonlocal_nu_nu_k_k_w_w);

          //SHOW::execute(MOMS_imag.G0_r_t);

          MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(MOMS_imag.G0_r_w, MOMS_imag.G0_k_w);
          MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(MOMS_real.G0_r_w, MOMS_real.G0_k_w);
          MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(MOMS_imag.G0_r_t, MOMS_imag.G0_k_t);

          //SHOW::execute(MOMS_imag.G0_r_w);
          //SHOW::execute(MOMS_imag.G0_k_w);
        }

      //assert(false);
    }

    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_real_space_Greens_functions(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
                                                                                                       FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
                                                                                                       FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                                                                                       FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> >& /*G_nu_nu_r_r_w_w*/,
                                                                                                       FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> >& /*G_nu_nu_k_k_w_w*/)
    {
      if(concurrency.id()==0)
        std::cout << "\n\t" << __FUNCTION__ << std::endl;

      G_r_w_im = 0;
      G_r_w_re = 0;
      G_r_t = 0;

      int start = clock();

      {
        int n_threads = 1;

        std::vector<sp_Greens_function_data_type> data_vec(n_threads);

        for(int l=0; l<n_threads; l++)
          data_vec[l].initialize(parameters);

        //         compute_Greens_functions_ac_slow(data_vec);//, G_r_w_im, G_r_w_re, G_r_t);
        //         compute_Greens_functions_ca_slow(data_vec);//, G_r_w_im, G_r_w_re, G_r_t);

        compute_Greens_functions_slow(data_vec);

        for(int l=0; l<n_threads; l++){
          data_vec[l].sum_to(G_r_w_im, G_r_w_re, G_r_t);
          //           data_vec[l].sum_to(G_nu_nu_r_r_w_w, G_nu_nu_k_k_w_w);
        }
      }

      {
        for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
          for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
            for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
              for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                G_r_t(nu_i, nu_j, r_i, t_i-t::dmn_size()/2) = -G_r_t(nu_i, nu_j, r_i, t_i);
      }

      int end = clock();

      if(concurrency.id()==0)
        std::cout << "\t" << __FUNCTION__ << " total time : " << double(end-start)/double(CLOCKS_PER_SEC) << std::endl;
    }

    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::renormalize_real_space_Greens_functions(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w_im,
                                                                                                           FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_re,
                                                                                                           FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                                                                                           FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, r_dmn, r_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_r_r_w_w,
                                                                                                           FUNC_LIB::function<complex_type, dmn_6<nu_dmn, nu_dmn, k_dmn, k_dmn, w_VERTEX, w_VERTEX> >& G_nu_nu_k_k_w_w)
    {
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n\t" << __FUNCTION__ << std::endl;
      }

      std::vector< Hilbert_space_type >& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      double beta = parameters.get_beta();

      double Z = 0;
      for(int HS_i = 0; HS_i < Hilbert_spaces.size(); ++HS_i)
        for(int n = 0; n < Hilbert_spaces[HS_i].size(); ++n)
          Z += std::exp(-beta*eigen_energies(HS_i)[n]);

      double factor = 1./Z;

      if (concurrency.id() == concurrency.first()) {
        std::cout << "\tfactor = " << factor << std::endl;
      }

      G_r_w_im *= factor;
      G_r_w_re *= factor;
      G_r_t    *= -factor;

      G_nu_nu_r_r_w_w *= factor;
      G_nu_nu_k_k_w_w *= factor;
    }

    /*
      template<typename parameter_type, typename ed_options>
      template<typename MOMS_w_imag_type, typename MOMS_w_real_type>
      void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_all_sp_functions(MOMS_w_imag_type& MOMS_imag,
      MOMS_w_real_type& MOMS_real,
      bool interacting)
      {
      if(interacting){
      compute_Greens_functions_main(MOMS_imag.G_r_w , MOMS_real.G_r_w , MOMS_imag.G_r_t,
      MOMS_imag.G_k_w , MOMS_real.G_k_w , MOMS_imag.G_k_t);
      }
      else{
      compute_Greens_functions_main(MOMS_imag.G0_r_w , MOMS_real.G0_r_w , MOMS_imag.G0_r_t,
      MOMS_imag.G0_k_w , MOMS_real.G0_k_w , MOMS_imag.G0_k_t);
      }
      }
    */

    /*
      template<typename parameter_type, typename ed_options>
      void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_Greens_functions_main(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
      FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w     > >& G_k_w,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_REAL> >& G_k_w_real,
      FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, k_dmn, t     > >& G_k_t)
      {
      if(concurrency.id()==0)
      std::cout << "\n\t" << __FUNCTION__ << std::endl;

      G_r_w = 0;
      G_k_w = 0;

      G_r_t = 0;
      G_k_t = 0;

      int start = clock();

      if(false)
      {
      compute_Greens_functions_ac_translation(G_r_w, G_r_w_real, G_r_t);
      compute_Greens_functions_ca_translation(G_r_w, G_r_w_real, G_r_t);
      }
      else
      {
      int n_threads = 1;

      std::vector<sp_Greens_function_data_type> data_vec(n_threads);

      for(int l=0; l<n_threads; l++)
      data_vec[l].initialize(parameters);

      compute_Greens_functions_ac_slow(data_vec, G_r_w, G_r_w_real, G_r_t);
      compute_Greens_functions_ca_slow(data_vec, G_r_w, G_r_w_real, G_r_t);
      }

      concurrency.sum(G_r_t);
      concurrency.sum(G_r_w);
      concurrency.sum(G_r_w_real);

      int end = clock();

      if(concurrency.id()==0)
      std::cout << "\t total time : " << double(end-start)/double(CLOCKS_PER_SEC) << "\n\n";

      {
      for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
      for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
      for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
      for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
      G_r_t(nu_i, nu_j, r_i, t_i-t::dmn_size()/2) = -G_r_t(nu_i, nu_j, r_i, t_i);
      }

      std::vector< Hilbert_space_type >& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      {
      double beta = parameters.get_beta();

      double Z = 0;
      for(int HS_i = 0; HS_i < Hilbert_spaces.size(); ++HS_i)
      for(int n = 0; n < Hilbert_spaces[HS_i].size(); ++n)
      Z += std::exp(-beta*eigen_energies(HS_i)[n]);

      //double factor = 1./(Z*r_dmn::dmn_size());
      double factor = 1./Z;

      std::cout << factor << std::endl;

      G_r_w      *= factor;
      G_r_w_real *= factor;
      G_r_t      *= -factor;

      //       FT<r_dmn, k_dmn>::execute(G_r_w     , G_k_w);
      //       FT<r_dmn, k_dmn>::execute(G_r_w_real, G_k_w_real);
      //       FT<r_dmn, k_dmn>::execute(G_r_t     , G_k_t);
      MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(G_r_w     , G_k_w);
      MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(G_r_w_real, G_k_w_real);
      MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(G_r_t     , G_k_t);

      SHOW::execute(G_r_t);
      //         SHOW::execute(G_k_t);

      //         SHOW::execute(G_r_w);
      SHOW::execute(G_k_w);
      }

      //assert(false);
      }
    */

    /*
      template<typename parameter_type, typename ed_options>
      void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_Greens_functions_ac_translation(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
      FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
      {
      std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

      std::vector< Hilbert_space_type >& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      scalar_type beta = parameters.get_beta();

      double zero=0, nonzero=0;

      #ifdef USE_OPENMP_DIRECTIVES
      //omp_set_num_threads(4);
      #pragma omp parallel for
      #endif
      for(int ind=0; ind<nu_nu_r_dmn_type::dmn_size(); ind++)
      {
      sp_Greens_function_data_type sp_Greens_function_data;//(parameters);

      sp_Greens_function_data.set_indices(ind);

      {
      stringstream ss;
      ss << "\t bsr_i : " << sp_Greens_function_data.bsr_i
      << "\t bsr_j : " << sp_Greens_function_data.bsr_j
      << "\t" << print_time()
      << "\n";
      std::cout << ss.str();
      }

      for(int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0)
      {
      for(int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1)
      {
      bool done = false;

      if(annihilation_set_all(HS_0, HS_1, sp_Greens_function_data.bsr_i) != -1 &&
      creation_set_all    (HS_1, HS_0, sp_Greens_function_data.bsr_j) != -1)
      {
      for(int n0 = 0; n0 < Hilbert_spaces[HS_0].size(); ++n0)
      {
      scalar_type E_n0 = eigen_energies(HS_0)[n0];
      scalar_type w_e  = std::exp(-beta*E_n0);

      if(w_e > CUT_OFF)
      {
      for(int n1 = 0; n1 < Hilbert_spaces[HS_1].size(); ++n1)
      {
      if(not done)
      {
      overlap.compute_annihilation_matrix_fast(HS_0, HS_1, sp_Greens_function_data.bsr_i, sp_Greens_function_data.annihilation_bsr_i, sp_Greens_function_data.tmp);
      overlap.compute_creation_matrix_fast    (HS_1, HS_0, sp_Greens_function_data.bsr_j, sp_Greens_function_data.creation_bsr_j    , sp_Greens_function_data.tmp);

      done = true;
      }

      scalar_type E_n1 = eigen_energies(HS_1)[n1];

      complex_type factor = w_e * sp_Greens_function_data.annihilation_bsr_i(n0,n1) * sp_Greens_function_data.creation_bsr_j(n1,n0);

      if(abs(factor) > CUT_OFF)
      {
      nonzero += 1;

      compute_sp_Greens_function(sp_Greens_function_data.nu_i, sp_Greens_function_data.nu_j, sp_Greens_function_data.delta_r,
      E_n0, E_n1, factor,
      G_r_w, G_r_w_real, G_r_t);
      }
      else
      {
      zero += 1;
      }
      }
      }
      }
      }
      }
      }
      }

      std::cout << "\n\t nonzero/total : " << nonzero/(zero+nonzero) << "\n";
      }
    */

    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_Greens_functions_ac_slow(std::vector<sp_Greens_function_data_type>& data_vec)//,
    //                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
    //                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
    //                                                                                                     FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

      std::vector< Hilbert_space_type >& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      scalar_type beta = parameters.get_beta();

      double zero=0, nonzero=0;

      for(int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0)
        {
          for(int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1)
            {
              for(int ind=0; ind<nu_nu_r_dmn_type::dmn_size(); ind++)
                {
                  sp_Greens_function_data_type& data = data_vec[0];

                  data.set_indices(ind);

                  bool done = false;

                  if(annihilation_set_all(HS_0, HS_1, data.bsr_i) != -1 &&
                     creation_set_all    (HS_1, HS_0, data.bsr_j) != -1)
                    {
                      for(int n0 = 0; n0 < Hilbert_spaces[HS_0].size(); ++n0)
                        {
                          scalar_type E_n0 = eigen_energies(HS_0)[n0];
                          scalar_type w_e  = std::exp(-beta*E_n0);

                          if(w_e > CUT_OFF)
                            {
                              for(int n1 = 0; n1 < Hilbert_spaces[HS_1].size(); ++n1)
                                {
                                  if(not done)
                                    {
                                      overlap.compute_annihilation_matrix_fast(HS_0, HS_1, data.bsr_i, data.annihilation_bsr_i, data.tmp);
                                      overlap.compute_creation_matrix_fast    (HS_1, HS_0, data.bsr_j, data.creation_bsr_j    , data.tmp);

                                      done = true;
                                    }

                                  scalar_type E_n1 = eigen_energies(HS_1)[n1];

                                  complex_type factor = w_e * data.annihilation_bsr_i(n0,n1) * data.creation_bsr_j(n1,n0);

                                  if(abs(factor) > CUT_OFF)
                                    {
                                      nonzero += 1;

                                      //                                       compute_sp_Greens_function(data.nu_i, data.nu_j, data.delta_r,
                                      //                                                                  E_n0, E_n1, factor,
                                      //                                                                  G_r_w, G_r_w_real, G_r_t);

                                      compute_sp_Greens_function(ind, E_n0, E_n1, factor, data);
                                    }
                                  else
                                    {
                                      zero += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

      std::cout << "\n\t nonzero/total : " << nonzero/(zero+nonzero) << "\n";
    }

    /*
      template<typename parameter_type, typename ed_options>
      void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_Greens_functions_ca_translation(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
      FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
      {
      std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

      std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      scalar_type beta = parameters.get_beta();

      double zero=0, nonzero=0;

      #ifdef USE_OPENMP_DIRECTIVES
      //omp_set_num_threads(4);
      #pragma omp parallel for
      #endif
      for(int ind=0; ind<nu_nu_r_dmn_type::dmn_size(); ind++)
      {
      sp_Greens_function_data_type sp_Greens_function_data;//(parameters);

      sp_Greens_function_data.set_indices(ind);

      {
      stringstream ss;
      ss << "\t bsr_i : " << sp_Greens_function_data.bsr_i
      << "\t bsr_j : " << sp_Greens_function_data.bsr_j
      << "\t" << print_time()
      << "\n";
      std::cout << ss.str();
      }

      for(int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0)
      {
      for(int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1)
      {
      bool done = false;

      if(annihilation_set_all(HS_1, HS_0, sp_Greens_function_data.bsr_i) != -1 &&
      creation_set_all    (HS_0, HS_1, sp_Greens_function_data.bsr_j) != -1)
      {
      for(int n0 = 0; n0 < Hilbert_spaces[HS_0].size(); ++n0)
      {
      scalar_type E_n0 = eigen_energies(HS_0)[n0];
      scalar_type w_e  = std::exp(-beta*E_n0);

      if(w_e > CUT_OFF)
      {
      for(int n1 = 0; n1 < Hilbert_spaces[HS_1].size(); ++n1)
      {
      if(not done)
      {
      overlap.compute_annihilation_matrix_fast(HS_1, HS_0, sp_Greens_function_data.bsr_i, sp_Greens_function_data.annihilation_bsr_i, sp_Greens_function_data.tmp);
      overlap.compute_creation_matrix_fast    (HS_0, HS_1, sp_Greens_function_data.bsr_j, sp_Greens_function_data.creation_bsr_j    , sp_Greens_function_data.tmp);

      done = true;
      }

      scalar_type E_n1 = eigen_energies(HS_1)[n1];

      complex_type factor = w_e * sp_Greens_function_data.creation_bsr_j(n0,n1) * sp_Greens_function_data.annihilation_bsr_i(n1,n0);

      if(abs(factor) > CUT_OFF)
      {
      nonzero += 1;

      compute_sp_Greens_function(sp_Greens_function_data.nu_i, sp_Greens_function_data.nu_j, sp_Greens_function_data.delta_r,
      E_n1, E_n0, factor,
      G_r_w, G_r_w_real, G_r_t);
      }
      else
      {
      zero += 1;
      }
      }
      }
      }
      }
      }
      }
      }

      std::cout << "\n\t nonzero/total : " << nonzero/(zero+nonzero) << "\n";
      }
    */

    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_Greens_functions_ca_slow(std::vector<sp_Greens_function_data_type>& data_vec)// ,
    //                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
    //                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
    //                                                                                                     FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      if(concurrency.id()==0)
        std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

      std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      scalar_type beta = parameters.get_beta();

      double zero=0, nonzero=0;

      for(int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0)
        {
          for(int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1)
            {
              for(int ind=0; ind<nu_nu_r_dmn_type::dmn_size(); ind++)
                {
                  sp_Greens_function_data_type& data = data_vec[0];

                  data.set_indices(ind);

                  bool done = false;

                  if(annihilation_set_all(HS_1, HS_0, data.bsr_i) != -1 &&
                     creation_set_all    (HS_0, HS_1, data.bsr_j) != -1)
                    {
                      for(int n0 = 0; n0 < Hilbert_spaces[HS_0].size(); ++n0)
                        {
                          scalar_type E_n0 = eigen_energies(HS_0)[n0];
                          scalar_type w_e  = std::exp(-beta*E_n0);

                          if(w_e > CUT_OFF)
                            {
                              for(int n1 = 0; n1 < Hilbert_spaces[HS_1].size(); ++n1)
                                {
                                  if(not done)
                                    {
                                      overlap.compute_annihilation_matrix_fast(HS_1, HS_0, data.bsr_i, data.annihilation_bsr_i, data.tmp);
                                      overlap.compute_creation_matrix_fast    (HS_0, HS_1, data.bsr_j, data.creation_bsr_j    , data.tmp);

                                      done = true;
                                    }

                                  scalar_type E_n1 = eigen_energies(HS_1)[n1];

                                  complex_type factor = w_e * data.creation_bsr_j(n0,n1) * data.annihilation_bsr_i(n1,n0);

                                  if(abs(factor) > CUT_OFF)
                                    {
                                      nonzero += 1;

                                      //                                       compute_sp_Greens_function(data.nu_i, data.nu_j, data.delta_r,
                                      //                                                                  E_n1, E_n0, factor,
                                      //                                                                  G_r_w, G_r_w_real, G_r_t);

                                      compute_sp_Greens_function(ind, E_n1, E_n0, factor, data);

                                    }
                                  else
                                    {
                                      zero += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

      std::cout << "\n\t nonzero/total : " << nonzero/(zero+nonzero) << "\n";
    }

    /*
      template<typename parameter_type, typename ed_options>
      void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_sp_Greens_function(int nu_i, int nu_j, int delta_r,
      scalar_type E_0, scalar_type E_1, complex_type factor,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
      FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
      {
      scalar_type beta    = parameters.get_beta();
      scalar_type off_set = parameters.get_real_frequencies_off_set();

      scalar_type  ONE(1);
      complex_type I(0,1);

      scalar_type  factor_re = real(factor);

      for(int t_i = t::dmn_size()/2; t_i < t::dmn_size(); ++t_i)
      {
      scalar_type tau   = t::get_elements()[t_i];
      scalar_type G_tau = 1./(std::exp((beta-tau)*(E_0 - E_1)) + std::exp(-tau*(E_0 - E_1)));

      //G_r_t(nu_i, nu_j, delta_r, t_i) += w_e * G_tau * real(annihilation_bsr_i(n0,n1)) * real(creation_bsr_j(n1,n0));
      G_r_t(nu_i, nu_j, delta_r, t_i) += factor_re * G_tau;// * real(annihilation_bsr_i(n0,n1)) * real(creation_bsr_j(n1,n0));
      }

      for(int w_i=0; w_i<w::dmn_size(); ++w_i)
      {
      scalar_type w_im = w::get_elements()[w_i];

      complex_type iw   = I*w_im;
      complex_type G_w  = ONE/(iw + E_0 - E_1);

      //G_r_w(nu_i, nu_j, delta_r, w_i) += w_e * G_w * annihilation_bsr_i(n0,n1) * creation_bsr_j(n1,n0);
      G_r_w(nu_i, nu_j, delta_r, w_i) += factor * G_w;// * annihilation_bsr_i(n0,n1) * creation_bsr_j(n1,n0);
      }

      for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i)
      {
      scalar_type w_re = w_REAL::get_elements()[w_i];

      complex_type iw  = w_re + I*off_set;
      complex_type G_w = ONE/(iw + E_0 - E_1);

      //G_r_w_real(nu_i, nu_j, delta_r, w_i) += w_e * G_w * annihilation_bsr_i(n0,n1) * creation_bsr_j(n1,n0);
      G_r_w_real(nu_i, nu_j, delta_r, w_i) += factor * G_w;// * annihilation_bsr_i(n0,n1) * creation_bsr_j(n1,n0);
      }
      }
    */

    /*!
     *   G(\tau) = e^(\tau*(E_0-E_1))/(e^(\beta*(E_0-E_1))+1)
     *   G(\omega) = 1/(i*\omega + (E_0-E_1))
     *   G(\tau) = 1/\beta \sum_{m=-\infty}^{\infty} G(\omega_m) e^{i*\omega_m*\tau}
     */
    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_sp_Greens_function(int nu_i_nu_j_delta_r,
                                                                                              scalar_type  E_0,
                                                                                              scalar_type  E_1,
                                                                                              complex_type factor,
                                                                                              sp_Greens_function_data_type& data)
    {
      scalar_type beta    = parameters.get_beta();

      scalar_type ONE(1);

      {
        for(int t_i = t::dmn_size()/2; t_i < t::dmn_size(); ++t_i)
          data.G_tau(t_i) = 1./(std::exp((beta-data.tau(t_i))*(E_0 - E_1)) + std::exp(-data.tau(t_i)*(E_0 - E_1)));

        for(int w_i=0; w_i<w::dmn_size(); ++w_i)
          data.G_w_im(w_i) = ONE/(data.w_im(w_i) + E_0 - E_1);

        for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i)
          data.G_w_re(w_i) = ONE/(data.w_re(w_i) + E_0 - E_1);
      }

      {
        scalar_type  factor_re = real(factor);

        for(int t_i = t::dmn_size()/2; t_i < t::dmn_size(); ++t_i)
          data.G_tau__nu_nu_r(t_i, nu_i_nu_j_delta_r) += factor_re * data.G_tau(t_i);

        for(int w_i=0; w_i<w::dmn_size(); ++w_i)
          data.G_w_im__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_im(w_i);

        for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i)
          data.G_w_re__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_re(w_i);
      }
    }

    template<typename parameter_type, typename ed_options>
    int fermionic_sp_Greens_function<parameter_type, ed_options>::has_nonzero_overlap(int HS_i, int HS_j,
                                                                                      bool is_creation, int bsr_ind)
    {
      if(is_creation)
        return creation_set_all(HS_i, HS_j, bsr_ind);
      else
        return annihilation_set_all(HS_i, HS_j, bsr_ind);
    }

    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::get_nonzero_overlap(int HS_i, int HS_j,
                                                                                       bool is_creation, int bsr_ind,
                                                                                       matrix_type& matrix,
                                                                                       matrix_type& tmp)
    {
      if(is_creation)
        overlap.compute_creation_matrix_fast    (HS_i, HS_j, bsr_ind, matrix, tmp);
      else
        overlap.compute_annihilation_matrix_fast(HS_i, HS_j, bsr_ind, matrix, tmp);
    }

    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_sp_permutations(int bsr_0, int bsr_1,
                                                                                           std::vector<std::vector<c_operator> >& sp_perms)
    {
      sp_perms.resize(0);

      std::vector<c_operator> c_operators(2);

      {
        c_operators[0].index = 0;
        c_operators[1].index = 1;

        c_operators[0].bsr_ind = bsr_0;
        c_operators[1].bsr_ind = bsr_1;

        c_operators[0].creation = false;
        c_operators[1].creation = true;
      }

      sp_perms.push_back(c_operators);

      //       {
      //         std::vector<int> indices(2);
      //         for(int l=0; l<2; l++)
      //           indices[l] = l;

      //         do
      //           {
      //             std::vector<c_operator> tmp;

      //             for(int l=0; l<2; l++)
      //               tmp.push_back(c_operators[indices[l]]);

      //             sp_perms.push_back(tmp);
      //           }
      //         while(std::next_permutation(indices.begin(), indices.end()));
      //       }
    }

    /*!
     *   G(\tau)   = e^{(E_i-E_j)\tau}
     *   G(\omega) =
     *   G(\tau)   =
     */
    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_sp_Greens_function_slow(int nu_i_nu_j_delta_r,
                                                                                                   scalar_type  E_i,
                                                                                                   scalar_type  E_j,
                                                                                                   complex_type factor,
                                                                                                   sp_Greens_function_data_type& data)
    {
      scalar_type beta = parameters.get_beta();

      complex_type e_min_beta_E_i = std::exp(-beta*E_i);
      complex_type e_min_beta_E_j = std::exp(-beta*E_j);

      if(abs(e_min_beta_E_i)>CUT_OFF or
         abs(e_min_beta_E_j)>CUT_OFF)
        {
          scalar_type ONE(1);

          {
            data.G_tau = 0;
            for(int t_i = t::dmn_size()/2; t_i < t::dmn_size(); ++t_i){
              data.G_tau(t_i) = std::exp(-E_i*(beta-data.tau(t_i)))*std::exp(-E_j*(data.tau(t_i)-0.));
            }

            data.G_w_im = 0;
            for(int w_i=0; w_i<w::dmn_size(); ++w_i){
              complex_type tmp = ONE/(data.w_im(w_i) + E_i - E_j);
              data.G_w_im(w_i) += std::exp(-beta*E_i)*tmp;
              data.G_w_im(w_i) += std::exp(-beta*E_j)*tmp;
            }

            data.G_w_re = 0;
            for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i){
              complex_type tmp = ONE/(data.w_re(w_i) + E_i - E_j);
              data.G_w_re(w_i) += std::exp(-beta*E_i)*tmp;
              data.G_w_re(w_i) += std::exp(-beta*E_j)*tmp;
            }

          }

          {
            scalar_type  factor_re = real(factor);

            for(int t_i = t::dmn_size()/2; t_i < t::dmn_size(); ++t_i)
              data.G_tau__nu_nu_r(t_i, nu_i_nu_j_delta_r) += factor_re * data.G_tau(t_i);

            for(int w_i=0; w_i<w::dmn_size(); ++w_i)
              data.G_w_im__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_im(w_i);

            for(int w_i=0; w_i<w_REAL::dmn_size(); ++w_i)
              data.G_w_re__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_re(w_i);
          }
        }
    }


    /*!
     *   G_{\nu, \mu}(R, \omega) = \frac{1}{Z} \sum_{i,j} e^{-\beta E_i} \phi(E_i, E_j, \omega) * \langle i | c_{\nu, \vec{R}} | j \rangle * \langle j | c^{\dagger}_{\mu, \vec{0}} | j \rangle
     *   \phi(E_i, E_j, \omega)  = \int_0^{\beta} d\tau e^{((E_i-E_j)+i \omega)\tau}
     *                           = \frac{e^{((E_i-E_j)+i \omega)\beta}-1}{((E_i-E_j)+i \omega)}
     */
    template<typename parameter_type, typename ed_options>
    void fermionic_sp_Greens_function<parameter_type, ed_options>::compute_Greens_functions_slow(std::vector<sp_Greens_function_data_type>& data_vec)
    {
      if(concurrency.id()==0)
        std::cout << "\t" << __FUNCTION__ << std::endl;;

      int origin = r_dmn::parameter_type::origin_index();

      //scalar_type beta = parameters.get_beta();

      std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

      for(int HS_0=0; HS_0<Hilbert_spaces.size(); ++HS_0){
        for(int HS_1=0; HS_1<Hilbert_spaces.size(); ++HS_1){

          for(int nu_0=0; nu_0<nu::dmn_size(); nu_0++){
            for(int nu_1=0; nu_1<nu::dmn_size(); nu_1++){

              for(int r_0=0; r_0<r_dmn::dmn_size(); r_0++){

                int r_1 = origin;
                //for(int r_1=0; r_1<r_dmn::dmn_size(); r_1++){

                sp_Greens_function_data_type& data = data_vec[0];

                int nu_0_nu_1_r_0 = data.nu_nu_r_dmn(nu_0, nu_1, r_0);

                int bsr_0 = data.nu_r_dmn(nu_0, r_0);
                int bsr_1 = data.nu_r_dmn(nu_1, r_1);

                std::vector<std::vector<c_operator> > sp_perms;

                compute_sp_permutations(bsr_0, bsr_1, sp_perms);

                for(int prm_ind=0; prm_ind<sp_perms.size(); prm_ind++){

                  scalar_type sign = 1;//-((prm_ind%2)-0.5)*2.0;

                  std::vector<c_operator>& operators = sp_perms[prm_ind];

                  if(has_nonzero_overlap(HS_0, HS_1, operators[0].creation, operators[0].bsr_ind) != -1 &&
                     has_nonzero_overlap(HS_1, HS_0, operators[1].creation, operators[1].bsr_ind) != -1){

                    bool done = false;

                    for(int l_0=0; l_0<Hilbert_spaces[HS_0].size(); ++l_0){

                      scalar_type E_0 = eigen_energies(HS_0)[l_0];

                      for(int l_1=0; l_1<Hilbert_spaces[HS_1].size(); ++l_1){

                        scalar_type E_1 = eigen_energies(HS_1)[l_1];

                        if(not done)
                          {
                            get_nonzero_overlap(HS_0, HS_1, operators[0].creation, operators[0].bsr_ind, data.overlap_0, data.tmp);
                            get_nonzero_overlap(HS_1, HS_0, operators[1].creation, operators[1].bsr_ind, data.overlap_1, data.tmp);

                            done = true;
                          }

                        complex_type factor = sign;//*std::exp(-beta*E_0);

                        factor *= data.overlap_0(l_0, l_1);
                        factor *= data.overlap_1(l_1, l_0);

                        if(abs(factor)>CUT_OFF){

                          compute_sp_Greens_function_slow(nu_0_nu_1_r_0, E_0, E_1, factor, data);

                          /*
                            int w[2];
                            int k[2];

                            for(int w2=0; w2<w_VERTEX::dmn_size(); w2++){
                            for(int w1=0; w1<w_VERTEX::dmn_size(); w1++){

                            if(abs(sign-1)>1.e-6){
                            w[operators[0].index] = w1;//w_VERTEX::dmn_size()-1-
                            w[operators[1].index] = w2;
                            }
                            else{
                            w[operators[0].index] = w_VERTEX::dmn_size()-1-w1;//w_VERTEX::dmn_size()-1-
                            w[operators[1].index] = w_VERTEX::dmn_size()-1-w2;
                            }

                            complex_type e_iwt = phi_e_iwt(E_0, E_1, w[0], w[1]);

                            data.G2_nonlocal_nu_nu_r_r_w_w(nu_0, nu_1, r_0, r_1, w1, w2) += sign*factor*e_iwt;

                            for(int k1=0; k1<k_dmn::dmn_size(); k1++){
                            for(int k2=0; k2<k_dmn::dmn_size(); k2++){

                            //                                     k[operators[0].index] = k_dmn::parameter_type::subtract(k1, origin);
                            //                                     k[operators[1].index] = k2;//k_dmn::parameter_type::subtract(k2, origin)

                            k[0] = k1;
                            k[1] = k_dmn::parameter_type::subtract(k2, origin);

                            complex_type e_ikr = phi_e_ikr(r_0, r_1, k[0], k[1]);

                            data.G2_nonlocal_nu_nu_k_k_w_w(nu_0, nu_1, k1, k2, w1, w2) += sign*factor*e_iwt*e_ikr;
                            }
                            }
                            }
                            }
                          */
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      //   }

      {
        //         sp_Greens_function_data_type& data = data_vec[0];

        //         int N = 8;

        //         std::cout << "\n\t G0_nonlocal_r_w\n\n";
        //         for(int w_i=0; w_i<w_VERTEX::dmn_size(); w_i++){
        //           std::cout << w_VERTEX::get_elements()[w_i] << "\t";
        //           for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
        //             std::cout << real(data.G2_nonlocal_nu_nu_r_r_w_w(0,0,k_i,0,w_i,w_i)) << "\t" << imag(data.G2_nonlocal_nu_nu_r_r_w_w(0,0,k_i,0,w_i,w_i)) << "\t";
        //           std::cout << "\n";
        //         }
        //         std::cout << "\n";

        //         std::cout << "\n\t G0_r_w\n\n";
        //         for(int w_i=w::dmn_size()/2-N; w_i<w::dmn_size()/2+N; w_i++){
        //           std::cout << w::get_elements()[w_i] << "\t";
        //           for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
        //             std::cout << real(data.G_w_im__nu_nu_r(w_i, 0,0, 0,0, k_i)) << "\t" << imag(data.G_w_im__nu_nu_r(w_i, 0,0, 0,0, k_i)) << "\t";
        //           std::cout << "\n";
        //         }
        //         std::cout << "\n";

        //         std::cout << "\n\t G0_nonlocal_k_w\n\n";
        //         for(int w_i=0; w_i<w_VERTEX::dmn_size(); w_i++){
        //           std::cout << w_VERTEX::get_elements()[w_i] << "\t";
        //           for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
        //             std::cout << real(data.G2_nonlocal_nu_nu_k_k_w_w(0,0,k_i,k_i,w_i,w_i)) << "\t" << imag(data.G2_nonlocal_nu_nu_k_k_w_w(0,0,k_i,k_i,w_i,w_i)) << "\t";
        //           std::cout << "\n";
        //         }
        //         std::cout << "\n";

        //         std::cout << "\n\t G0_k_w\n\n";
        //         for(int w_i=w::dmn_size()/2-N; w_i<w::dmn_size()/2+N; w_i++){
        //           std::cout << w::get_elements()[w_i] << "\t";
        //           for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
        //             std::cout << real(data.G_w_im__nu_nu_k(w_i, 0, 0, k_i)) << "\t" << imag(data.G_w_im__nu_nu_k(w_i, 0, 0, k_i)) << "\t";
        //           std::cout << "\n";
        //         }
        //         std::cout << "\n";
      }

      //       assert(false);
    }

  }

}

#endif
