//-*-C++-*-

#ifndef FERMIONIC_SP_GREENS_FUNCTION_H
#define FERMIONIC_SP_GREENS_FUNCTION_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    struct psi_lhs_domain
    {
      typedef int            element_type;
      typedef psi_lhs_domain this_type;

      static int& get_size()
      {
        static int size = 0;
        return size;
      }

      static std::vector<int>& get_elements()
      {
        static std::vector<int> elements(get_size(), 0);
        return elements;
      }
    };

    struct psi_rhs_domain
    {
      typedef int            element_type;
      typedef psi_rhs_domain this_type;

      static int& get_size()
      {
        static int size = 0;
        return size;
      }

      static std::vector<int>& get_elements()
      {
        static std::vector<int> elements(get_size(), 0);
        return elements;
      }
    };

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    class fermionic_sp_Greens_function
    {
#include "type_definitions.h"

      typedef ED_type_definitions<parameter_type, b_dmn, s_dmn, r_dmn> ED_type_def;

      typedef typename ED_type_def::profiler_t       profiler_t;
      typedef typename ED_type_def::concurrency_type concurrency_type;

      typedef typename ED_type_def::scalar_type  scalar_type;
      typedef typename ED_type_def::complex_type complex_type;

      typedef typename ED_type_def::vector_type         vector_type;
      typedef typename ED_type_def::matrix_type         matrix_type;
      typedef typename ED_type_def::int_matrix_type int_matrix_type;

      typedef typename ED_type_def::nu_dmn nu_dmn;

      typedef typename ED_type_def::occ_dmn occ_dmn;
      typedef typename ED_type_def::mag_dmn mag_dmn;

      typedef typename ED_type_def::occ_mag_dmn occ_mag_dmn;

      typedef fermionic_Fock_space <parameter_type, b_dmn, s_dmn, r_dmn> fermionic_Fock_space_type;
      typedef fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn> fermionic_Hamiltonian_type;

    public:

      fermionic_sp_Greens_function(parameter_type&             parameters_ref,
                                   fermionic_Fock_space_type&  Fock_space_ref,
                                   fermionic_Hamiltonian_type& fermionic_Hamiltonian_ref);

      ~fermionic_sp_Greens_function();

      template<typename MOMS_w_imag_type, typename MOMS_w_real_type>
      void compute_all_sp_functions(MOMS_w_imag_type& MOMS_imag, MOMS_w_real_type& MOMS_real, bool interacting);

      template<typename k_dmn, typename w_dmn>
      void compute_S_k_w(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G_k_w,
                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G0_k_w,
                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& S_k_w);

      template<typename MOMS_type>
      void print_G_k_w(MOMS_type& MOMS, bool interacting);

    private:

      template<typename k_dmn>
      void compute_Greens_functions(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                    FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                    FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                    FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w     > >& G_k_w,
                                    FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_REAL> >& G_k_w_real,
                                    FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, k_dmn, t     > >& G_k_t);

      void compute_ca_Greens_function_st(int n_0, int Sz_0,
                                         int n_1, int Sz_1,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                         FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_ca_Greens_function_mt(int n_0, int Sz_0,
                                         int n_1, int Sz_1,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                         FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_ac_Greens_function_st(int n_0, int Sz_0,
                                         int n_1, int Sz_1,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                         FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_ac_Greens_function_mt(int n_0, int Sz_0,
                                         int n_1, int Sz_1,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                         FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                         FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t);

      void compute_overlap_serial_slow(int n_0, int Sz_0,
                                       int n_1, int Sz_1,
                                       std::vector<overlap_indices>& overlap_lhs,
                                       std::vector<overlap_indices>& overlap_rhs,
                                       int l0 , int l1);

      void compute_overlap_serial_fast(int n_0, int Sz_0,
                                       int n_1, int Sz_1,
                                       std::vector<overlap_indices>& overlap_lhs,
                                       std::vector<overlap_indices>& overlap_rhs,
                                       std::vector<int>& bp_lhs,
                                       std::vector<int>& bp_rhs,
                                       int l0 , int l1);

      void set_up_psi_domain(int n_0, int Sz_0,
                             int n_1, int Sz_1);

      bool check_overlap_ac(int n_0, int Sz_0, int n_1, int Sz_1, overlap_indices& overlap_i, overlap_indices& overlap_j);
      bool check_overlap_ca(int n_0, int Sz_0, int n_1, int Sz_1, overlap_indices& overlap_i, overlap_indices& overlap_j);

    private:

      parameter_type&   parameters;
      concurrency_type& concurrency;

      double CUT_OFF;

      fermionic_Fock_space_type&  Fock_space;
      fermionic_Hamiltonian_type& Hamiltonian;


      FUNC_LIB::function<int            , dmn_2<occ_dmn, mag_dmn> >& n_occupation_states;
      FUNC_LIB::function<int_matrix_type, dmn_2<occ_dmn, mag_dmn> >&   occupation_states;

      FUNC_LIB::function<std::vector<overlap_indices>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& creation_overlap_of_states;
      FUNC_LIB::function<std::vector<overlap_indices>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& annihilation_overlap_of_states;

      FUNC_LIB::function<std::vector<int>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& creation_overlap_break_points;
      FUNC_LIB::function<std::vector<int>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& annihilation_overlap_break_points;

      FUNC_LIB::function<vector_type, dmn_2<occ_dmn, mag_dmn> >& eigen_energies;
      FUNC_LIB::function<matrix_type, dmn_2<occ_dmn, mag_dmn> >& eigen_states;

      FUNC_LIB::function<int, dmn_2<r_dmn, r_dmn> > rj_minus_ri;

      FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >         overlap;
      FUNC_LIB::function<std::complex<double>, dmn_2<dmn_3<b_dmn, s_dmn, r_dmn>, dmn_3<b_dmn, s_dmn, r_dmn> > > overlap_r_r;

      LIN_ALG::vector<int         , LIN_ALG::CPU> V_lhs_index;
      LIN_ALG::vector<int         , LIN_ALG::CPU> V_rhs_index;

      LIN_ALG::vector<complex_type, LIN_ALG::CPU> V_lhs_value;
      LIN_ALG::vector<complex_type, LIN_ALG::CPU> V_rhs_value;

      LIN_ALG::vector<complex_type, LIN_ALG::CPU> psi_0_vec;
      LIN_ALG::vector<complex_type, LIN_ALG::CPU> psi_1_vec;
    };

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::fermionic_sp_Greens_function(parameter_type&             parameters_ref,
                                                                                                    fermionic_Fock_space_type&  Fock_space_ref,
                                                                                                    fermionic_Hamiltonian_type& Hamiltonian_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CUT_OFF(parameters.get_eigenvalue_cut_off()),//1.e-6),

      Fock_space (Fock_space_ref),
      Hamiltonian(Hamiltonian_ref),

      n_occupation_states(Fock_space.n_occupation_states),
      occupation_states  (Fock_space.  occupation_states),

      creation_overlap_of_states    (Fock_space.creation_overlap_of_states),
      annihilation_overlap_of_states(Fock_space.annihilation_overlap_of_states),

      creation_overlap_break_points    (Fock_space.creation_overlap_break_points),
      annihilation_overlap_break_points(Fock_space.annihilation_overlap_break_points),

      eigen_energies(Hamiltonian.get_eigen_energies()),
      eigen_states  (Hamiltonian.get_eigen_states()),

      rj_minus_ri("rj_minus_ri"),

      overlap    ("overlap"),
      overlap_r_r("overlap_r_r"),

      V_lhs_index("V_lhs_index"),
      V_rhs_index("V_rhs_index"),

      V_lhs_value("V_lhs_value"),
      V_rhs_value("V_rhs_value"),

      psi_0_vec("psi_0_vec"),
      psi_1_vec("psi_1_vec")
    {
      for(int ri=0; ri<r_dmn::dmn_size(); ri++)
        for(int rj=0; rj<r_dmn::dmn_size(); rj++)
          rj_minus_ri(ri, rj) = r_dmn::parameter_type::subtract(ri, rj);
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::~fermionic_sp_Greens_function()
    {}

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    template<typename k_dmn, typename w_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_S_k_w(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G_k_w,
                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& G0_k_w,
                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_dmn> >& S_k_w)
    {
      if(concurrency.id()==0)
        cout << "\n\t" << __FUNCTION__ << endl;

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

      delete [] G_inverted_matrix;
      delete [] G0_cluster_excluded_inverted_matrix;
      delete [] Sigma_matrix;

      if(concurrency.id()==0)
        {
          cout << "\n";
          for(int w_i=w_dmn::dmn_size()/2-16; w_i<w_dmn::dmn_size()/2+16; w_i++){
            cout << w_dmn::get_elements()[w_i] << "\t";
            for(int k_i=0; k_i<k_dmn::dmn_size(); k_i++)
              cout << real(S_k_w(0,0,0,0,k_i,w_i)) << "\t" << imag(S_k_w(0,0,0,0,k_i,w_i)) << "\t";
            cout << "\n";
          }
          cout << "\n";
        }

    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    template<typename MOMS_type>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::print_G_k_w(MOMS_type& MOMS, bool interacting)
    {
      cout << "\n\t" << __FUNCTION__ << endl;

      if(interacting)
        {
          cout << "\n";
          for(int w_i=w::dmn_size()/2-36; w_i<w::dmn_size()/2+36; w_i++){
            cout << w::get_elements()[w_i] << "\t";
            for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++)
              cout << real(MOMS.G_k_w(0,0,0,0,k_i,w_i)) << "\t" << imag(MOMS.G_k_w(0,0,0,0,k_i,w_i)) << "\t";
            cout << "\n";
          }
          cout << "\n";
        }
      else
        {
          cout << "\n";
          for(int w_i=w::dmn_size()/2-16; w_i<w::dmn_size()/2+16; w_i++){
            cout << w::get_elements()[w_i] << "\t";
            for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++)
              cout << real(MOMS.G0_r_w(0,0,0,0,k_i,w_i)) << "\t" << imag(MOMS.G0_r_w(0,0,0,0,k_i,w_i)) << "\t";
            cout << "\n";
          }
          cout << "\n";
        }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    template<typename MOMS_w_imag_type, typename MOMS_w_real_type>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_all_sp_functions(MOMS_w_imag_type& MOMS_imag,
                                                                                                     MOMS_w_real_type& MOMS_real,
                                                                                                     bool interacting)
    {
      if(interacting)
        {
          compute_Greens_functions(MOMS_imag.G_r_w , MOMS_real.G_r_w , MOMS_imag.G_r_t,
                                   MOMS_imag.G_k_w , MOMS_real.G_k_w , MOMS_imag.G_k_t);
        }
      else
        {
          compute_Greens_functions(MOMS_imag.G0_r_w , MOMS_real.G0_r_w , MOMS_imag.G0_r_t,
                                   MOMS_imag.G0_k_w , MOMS_real.G0_k_w , MOMS_imag.G0_k_t);
        }
    }

    /*! p 134
     *
     *   G(\nu, \mu, z) = \frac{1}{Z} \sum_{n_0, n_1} \frac{\langle n_0 | c_{\nu} | n_1 \rangle \langle n_1 | c^{\dagger}_{\nu} | n_0 \rangle}{z+(E_n-E_{n'})}
     *
     *   G(\tau, \epsilon) = \frac{1}{\beta} \sum_{m} \frac{1}{i \varpi_m + \epsilon} e^{i \varpi \tau}
     *                     = \frac{1}{\beta} \sum_{m} (\frac{1}{i \varpi_m + \epsilon} -\frac{1}{i \varpi_m }) e^{i \varpi \tau} + 0.5
     *                     = (1 - \frac{1}{e^{-\beta \epsilon}+1} ) * e^{(\beta-\tau) \epsilon}
     *                     = \frac{1}{e^{-(\beta-\tau) \epsilon} + e^{\tau \epsilon} }
     *
     */
    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    template<typename k_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_Greens_functions(FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                                                                     FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t,
                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w     > >& G_k_w,
                                                                                                     FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, k_dmn, w_REAL> >& G_k_w_real,
                                                                                                     FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, k_dmn, t     > >& G_k_t)
    {
      if(concurrency.id()==0)
        cout << "\n\t" << __FUNCTION__ << endl;

      G_r_w = 0;
      G_k_w = 0;

      G_r_t = 0;
      G_k_t = 0;

      int start = clock();

      for(int n_0=0; n_0<occ_dmn::dmn_size(); n_0++){
        for(int Sz_0=0; Sz_0<mag_dmn::dmn_size(); Sz_0++){
          for(int n_1=0; n_1<occ_dmn::dmn_size(); n_1++){
            for(int Sz_1=0; Sz_1<mag_dmn::dmn_size(); Sz_1++){

              if(annihilation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()>0 and
                 creation_overlap_of_states    (n_1, Sz_1, n_0, Sz_0).size()>0)
                {
                  set_up_psi_domain(n_0, Sz_0, n_1, Sz_1);

                  //compute_ac_Greens_function_st(n_0, Sz_0, n_1, Sz_1, G_r_w, G_r_w_real, G_r_t);
                  compute_ac_Greens_function_mt(n_0, Sz_0, n_1, Sz_1, G_r_w, G_r_w_real, G_r_t);
                }

              if(creation_overlap_of_states    (n_0, Sz_0, n_1, Sz_1).size()>0 and
                 annihilation_overlap_of_states(n_1, Sz_1, n_0, Sz_0).size()>0)
                {
                  set_up_psi_domain(n_0, Sz_0, n_1, Sz_1);

                  //compute_ca_Greens_function_st(n_0, Sz_0, n_1, Sz_1, G_r_w, G_r_w_real, G_r_t);
                  compute_ca_Greens_function_mt(n_0, Sz_0, n_1, Sz_1, G_r_w, G_r_w_real, G_r_t);
                }
            }
          }
        }
      }

      concurrency.sum(G_r_t);
      concurrency.sum(G_r_w);
      concurrency.sum(G_r_w_real);

      int end = clock();

      if(concurrency.id()==0)
        cout << "\t total time : " << double(end-start)/double(CLOCKS_PER_SEC) << "\n\n";

      {
        for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
          for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
            for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
              for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                G_r_t(nu_i, nu_j, r_i, t_i-t::dmn_size()/2) = -G_r_t(nu_i, nu_j, r_i, t_i);
      }

      {
        double beta = parameters.get_beta();

        double Z = 0;
        for(int i=0; i<occ_dmn::dmn_size(); i++)
          for(int j=0; j<mag_dmn::dmn_size(); j++)
            for(int n=0; n<n_occupation_states(i,j); n++)
              Z += std::exp(-beta*eigen_energies(i,j)[n]);

        double factor = 1./(Z*r_dmn::dmn_size());

        G_r_w      *= factor;
        G_r_w_real *= factor;
        G_r_t      *= -factor;

        //       FT<r_dmn, k_dmn>::execute(G_r_w     , G_k_w);
        //       FT<r_dmn, k_dmn>::execute(G_r_w_real, G_k_w_real);
        //       FT<r_dmn, k_dmn>::execute(G_r_t     , G_k_t);
        MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(G_r_w     , G_k_w);
        MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(G_r_w_real, G_k_w_real);
        MATH_ALGORITHMS::TRANSFORM<r_dmn, k_dmn>::execute(G_r_t     , G_k_t);

      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_ac_Greens_function_st(int n_0, int Sz_0,
                                                                                                          int n_1, int Sz_1,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                                                                          FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      std::complex<double> I(0,1);

      double beta = parameters.get_beta();
      double off_set = parameters.get_real_frequencies_off_set();

      int N_0 = n_occupation_states(n_0, Sz_0);
      int N_1 = n_occupation_states(n_1, Sz_1);

      for(int l0=0; l0<N_0; l0++){

        double E_n0 = eigen_energies(n_0, Sz_0)[l0];
        double w_e  = std::exp(-beta*E_n0);

        if(w_e > CUT_OFF)
          {
            for(int l1=0; l1<N_1; l1++)
              {
                double E_n1 = eigen_energies(n_1, Sz_1)[l1];

                {
                  std::vector<overlap_indices>& overlap_lhs = annihilation_overlap_of_states(n_0, Sz_0, n_1, Sz_1);
                  std::vector<overlap_indices>& overlap_rhs = creation_overlap_of_states    (n_1, Sz_1, n_0, Sz_0);

                  std::vector<int>& bp_lhs = annihilation_overlap_break_points(n_0, Sz_0, n_1, Sz_1);
                  std::vector<int>& bp_rhs = creation_overlap_break_points    (n_1, Sz_1, n_0, Sz_0);

                  compute_overlap_serial_fast(n_0, Sz_0, n_1, Sz_1, overlap_lhs, overlap_rhs, bp_lhs, bp_rhs, l0, l1);
                }

                for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
                  {
                    double tau = t::get_elements()[t_i];

                    double G_tau = 1./(std::exp((beta-tau)*(E_n0 - E_n1)) + std::exp(-tau*(E_n0 - E_n1)));

                    for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                      for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                        for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                          G_r_t(nu_i, nu_j, r_i, t_i) += w_e*G_tau*real(overlap(nu_i, nu_j, r_i));
                  }

                for(int w_i=0; w_i<w::dmn_size(); w_i++)
                  {
                    std::complex<double> iw = I*w::get_elements()[w_i];

                    std::complex<double> G_w = 1./(iw + E_n0 - E_n1);

                    for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                      for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                        for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                          G_r_w(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
                  }

                for(int w_i=0; w_i<w_REAL::dmn_size(); w_i++)
                  {
                    std::complex<double> iw = w_REAL::get_elements()[w_i]+I*off_set;

                    std::complex<double> G_w = 1./(iw + E_n0 - E_n1);

                    for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                      for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                        for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                          G_r_w_real(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
                  }
              }
          }
      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_ac_Greens_function_mt(int n_0, int Sz_0,
                                                                                                          int n_1, int Sz_1,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                                                                          FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      std::complex<double> I(0,1);

      double beta    = parameters.get_beta();
      double off_set = parameters.get_real_frequencies_off_set();

      dmn_2<dmn_0<psi_lhs_domain>, dmn_0<psi_rhs_domain> > dmn;
      thread_manager_sum<concurrency_type>                 sum_manager(concurrency);

      //do
      {
        std::pair<int, int> bounds = sum_manager.get_bounds(dmn);

        int* coor = new int[2];

        for(int l=bounds.first; l<bounds.second; l++)
          {
            dmn.linind_2_subind(l, coor);

            int l0 = psi_lhs_domain::get_elements()[coor[0]];
            int l1 = psi_rhs_domain::get_elements()[coor[1]];

            double E_n0 = eigen_energies(n_0, Sz_0)[l0];
            double E_n1 = eigen_energies(n_1, Sz_1)[l1];

            double w_e  = std::exp(-beta*E_n0);
            assert(w_e > CUT_OFF);

            {
              std::vector<overlap_indices>& overlap_lhs = annihilation_overlap_of_states(n_0, Sz_0, n_1, Sz_1);
              std::vector<overlap_indices>& overlap_rhs = creation_overlap_of_states    (n_1, Sz_1, n_0, Sz_0);

              std::vector<int>& bp_lhs = annihilation_overlap_break_points(n_0, Sz_0, n_1, Sz_1);
              std::vector<int>& bp_rhs = creation_overlap_break_points    (n_1, Sz_1, n_0, Sz_0);

              compute_overlap_serial_fast(n_0, Sz_0, n_1, Sz_1, overlap_lhs, overlap_rhs, bp_lhs, bp_rhs, l0, l1);
            }

            for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
              {
                double tau = t::get_elements()[t_i];

                double G_tau = 1./(std::exp((beta-tau)*(E_n0 - E_n1)) + std::exp(-tau*(E_n0 - E_n1)));

                for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                  for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                    for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                      G_r_t(nu_i, nu_j, r_i, t_i) += w_e*G_tau*real(overlap(nu_i, nu_j, r_i));
              }

            for(int w_i=0; w_i<w::dmn_size(); w_i++)
              {
                std::complex<double> iw = I*w::get_elements()[w_i];

                std::complex<double> G_w = 1./(iw + E_n0 - E_n1);

                for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                  for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                    for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                      G_r_w(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
              }

            for(int w_i=0; w_i<w_REAL::dmn_size(); w_i++)
              {
                std::complex<double> iw = w_REAL::get_elements()[w_i]+I*off_set;

                std::complex<double> G_w = 1./(iw + E_n0 - E_n1);

                for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                  for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                    for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                      G_r_w_real(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
              }
          }

        delete [] coor;
      }
      //     while(!sum_manager.sum_and_check(G_r_t) and
      //          !sum_manager.sum_and_check(G_r_w) and
      //          !sum_manager.sum_and_check(G_r_w_real));
    }


    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_ca_Greens_function_st(int n_0, int Sz_0,
                                                                                                          int n_1, int Sz_1,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                                                                          FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      double beta    = parameters.get_beta();
      double off_set = parameters.get_real_frequencies_off_set();

      int N_0 = n_occupation_states(n_0, Sz_0);
      int N_1 = n_occupation_states(n_1, Sz_1);

      std::complex<double> I(0,1);

      for(int l0=0; l0<N_0; l0++){

        double E_n0 = eigen_energies(n_0, Sz_0)[l0];
        double w_e   = std::exp(-beta*E_n0);

        if(w_e > CUT_OFF)
          {
            for(int l1=0; l1<N_1; l1++){

              double E_n1 = eigen_energies(n_1, Sz_1)[l1];

              {
                std::vector<overlap_indices>& overlap_lhs = creation_overlap_of_states    (n_0, Sz_0, n_1, Sz_1);
                std::vector<overlap_indices>& overlap_rhs = annihilation_overlap_of_states(n_1, Sz_1, n_0, Sz_0);

                std::vector<int>& bp_lhs = creation_overlap_break_points    (n_0, Sz_0, n_1, Sz_1);
                std::vector<int>& bp_rhs = annihilation_overlap_break_points(n_1, Sz_1, n_0, Sz_0);

                compute_overlap_serial_fast(n_0, Sz_0, n_1, Sz_1, overlap_lhs, overlap_rhs, bp_lhs, bp_rhs, l0, l1);
              }

              for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
                {
                  double tau = t::get_elements()[t_i];

                  double G_tau = 1./(std::exp((beta-tau)*(E_n1 - E_n0)) + std::exp(-tau*(E_n1 - E_n0)));

                  for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                    for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                      for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                        G_r_t(nu_i, nu_j, r_i, t_i) += w_e*G_tau*real(overlap(nu_i, nu_j, r_i));
                }

              for(int w_i=0; w_i<w::dmn_size(); w_i++)
                {
                  std::complex<double> iw = I*w::get_elements()[w_i];

                  std::complex<double> G_w = 1./(iw - E_n0 + E_n1);

                  for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                    for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                      for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                        G_r_w(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
                }

              for(int w_i=0; w_i<w_REAL::dmn_size(); w_i++)
                {
                  std::complex<double> iw = w_REAL::get_elements()[w_i]+I*off_set;

                  std::complex<double> G_w = 1./(iw - E_n0 + E_n1);

                  for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                    for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                      for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                        G_r_w_real(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
                }
            }
          }
      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_ca_Greens_function_mt(int n_0, int Sz_0,
                                                                                                          int n_1, int Sz_1,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w     > >& G_r_w,
                                                                                                          FUNC_LIB::function<std::complex<double>, dmn_4<nu_dmn, nu_dmn, r_dmn, w_REAL> >& G_r_w_real,
                                                                                                          FUNC_LIB::function<             double , dmn_4<nu_dmn, nu_dmn, r_dmn, t     > >& G_r_t)
    {
      std::complex<double> I(0,1);

      double beta = parameters.get_beta();
      double off_set = parameters.get_real_frequencies_off_set();

      dmn_2<dmn_0<psi_lhs_domain>, dmn_0<psi_rhs_domain> > dmn;
      thread_manager_sum<concurrency_type>                 sum_manager(concurrency);

      //do
      {
        std::pair<int, int> bounds = sum_manager.get_bounds(dmn);

        int* coor = new int[2];

        for(int l=bounds.first; l<bounds.second; l++)
          {
            dmn.linind_2_subind(l, coor);

            int l0 = psi_lhs_domain::get_elements()[coor[0]];
            int l1 = psi_rhs_domain::get_elements()[coor[1]];

            double E_n0 = eigen_energies(n_0, Sz_0)[l0];
            double E_n1 = eigen_energies(n_1, Sz_1)[l1];

            double w_e  = std::exp(-beta*E_n0);
            assert(w_e > CUT_OFF);

            {
              std::vector<overlap_indices>& overlap_lhs = creation_overlap_of_states    (n_0, Sz_0, n_1, Sz_1);
              std::vector<overlap_indices>& overlap_rhs = annihilation_overlap_of_states(n_1, Sz_1, n_0, Sz_0);

              std::vector<int>& bp_lhs = creation_overlap_break_points    (n_0, Sz_0, n_1, Sz_1);
              std::vector<int>& bp_rhs = annihilation_overlap_break_points(n_1, Sz_1, n_0, Sz_0);

              compute_overlap_serial_fast(n_0, Sz_0, n_1, Sz_1, overlap_lhs, overlap_rhs, bp_lhs, bp_rhs, l0, l1);
            }

            for(int t_i=t::dmn_size()/2; t_i<t::dmn_size(); t_i++)
              {
                double tau = t::get_elements()[t_i];

                double G_tau = 1./(std::exp((beta-tau)*(E_n1 - E_n0)) + std::exp(-tau*(E_n1 - E_n0)));

                for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                  for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                    for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                      G_r_t(nu_i, nu_j, r_i, t_i) += w_e*G_tau*real(overlap(nu_i, nu_j, r_i));
              }

            for(int w_i=0; w_i<w::dmn_size(); w_i++)
              {
                std::complex<double> iw = I*w::get_elements()[w_i];

                std::complex<double> G_w = 1./(iw - E_n0 + E_n1);

                for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                  for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                    for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                      G_r_w(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
              }

            for(int w_i=0; w_i<w_REAL::dmn_size(); w_i++)
              {
                std::complex<double> iw = w_REAL::get_elements()[w_i]+I*off_set;

                std::complex<double> G_w = 1./(iw - E_n0 + E_n1);

                for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++)
                  for(int nu_j=0; nu_j<2*b_dmn::dmn_size(); nu_j++)
                    for(int nu_i=0; nu_i<2*b_dmn::dmn_size(); nu_i++)
                      G_r_w_real(nu_i, nu_j, r_i, w_i) += w_e*G_w*overlap(nu_i, nu_j, r_i);
              }
          }

        delete [] coor;
      }
      //     while(!sum_manager.sum_and_check(G_r_t) and
      //          !sum_manager.sum_and_check(G_r_w) and
      //          !sum_manager.sum_and_check(G_r_w_real));
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_overlap_serial_slow(int n_0, int Sz_0,
                                                                                                        int n_1, int Sz_1,
                                                                                                        std::vector<overlap_indices>& overlap_lhs,
                                                                                                        std::vector<overlap_indices>& overlap_rhs,
                                                                                                        int l0 , int l1)
    {
      overlap     = 0.;
      overlap_r_r = 0.;

      int N_lhs = overlap_lhs.size();
      int N_rhs = overlap_rhs.size();

      matrix_type& psi_0 = eigen_states(n_0, Sz_0);
      matrix_type& psi_1 = eigen_states(n_1, Sz_1);

      for(size_t j=0; j<N_rhs; j++){

        overlap_indices& overlap_j = overlap_rhs[j];

        for(size_t i=0; i<N_lhs; i++){

          overlap_indices& overlap_i = overlap_lhs[i];

          scalar_type phase = overlap_i.sign*overlap_j.sign;

          complex_type c_psi_0_i_times_psi_0_j = conj_value(psi_0(overlap_i.lhs, l0))*           psi_0(overlap_j.rhs, l0);
          complex_type psi_1_i_times_c_psi_0_j =            psi_1(overlap_i.rhs, l1) *conj_value(psi_1(overlap_j.lhs, l1));

          overlap_r_r(overlap_i.index, overlap_j.index) += phase*c_psi_0_i_times_psi_0_j*psi_1_i_times_c_psi_0_j;
        }
      }

      for(int r_j=0; r_j<r_dmn::dmn_size(); r_j++){
        for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++){

          int delta_r = rj_minus_ri(r_i, r_j);

          for(int s_j=0; s_j<s_dmn::dmn_size(); s_j++)
            for(int b_j=0; b_j<b_dmn::dmn_size(); b_j++)
              for(int s_i=0; s_i<s_dmn::dmn_size(); s_i++)
                for(int b_i=0; b_i<b_dmn::dmn_size(); b_i++)
                  overlap(b_i, s_i, b_j, s_j, delta_r)
                    += overlap_r_r(b_i, s_i, r_i, b_j, s_j, r_j);
        }
      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::compute_overlap_serial_fast(int n_0, int Sz_0,
                                                                                                        int n_1, int Sz_1,
                                                                                                        std::vector<overlap_indices>& overlap_lhs,
                                                                                                        std::vector<overlap_indices>& overlap_rhs,
                                                                                                        std::vector<int>& bp_lhs,
                                                                                                        std::vector<int>& bp_rhs,
                                                                                                        int l0 , int l1)
    {
      overlap     = 0.;
      overlap_r_r = 0.;

      matrix_type& psi_0 = eigen_states(n_0, Sz_0);
      matrix_type& psi_1 = eigen_states(n_1, Sz_1);

      psi_0_vec.reserve(psi_0.get_number_of_rows());
      psi_1_vec.reserve(psi_1.get_number_of_rows());

      memcpy(&psi_0_vec[0], &psi_0(0, l0), sizeof(complex_type)*psi_0.get_number_of_rows());
      memcpy(&psi_1_vec[0], &psi_1(0, l1), sizeof(complex_type)*psi_1.get_number_of_rows());

      int Nc_lhs = bp_lhs.size()-1;

      V_lhs_index.reserve(Nc_lhs);
      V_lhs_value.reserve(Nc_lhs);

      {// compute V_lhs
        for(size_t bp_i=0; bp_i<bp_lhs.size()-1; bp_i++){

          int delta_i = bp_lhs[bp_i+1]-bp_lhs[bp_i];
          int index_i = overlap_lhs[bp_i*delta_i].index;

          V_lhs_value[bp_i] = 0.;
          V_lhs_index[bp_i] = index_i;

          for(int o_i=0; o_i<delta_i; o_i++){

            int i = o_i + bp_i*delta_i;

            overlap_indices& overlap_i = overlap_lhs[i];
            assert(index_i == overlap_i.index);

            //V_lhs_value[bp_i] += overlap_i.sign*conj_value(psi_0(overlap_i.lhs, l0))*psi_1(overlap_i.rhs, l1);
            V_lhs_value[bp_i] += overlap_i.sign*conj_value(psi_0_vec[overlap_i.lhs])*psi_1_vec[overlap_i.rhs];
          }
        }
      }

      int Nc_rhs = bp_rhs.size()-1;

      V_rhs_index.reserve(Nc_rhs);
      V_rhs_value.reserve(Nc_rhs);

      {// compute V_rhs
        for(size_t bp_j=0; bp_j<bp_rhs.size()-1; bp_j++){

          int delta_j = bp_rhs[bp_j+1]-bp_rhs[bp_j];
          int index_j = overlap_rhs[bp_j*delta_j].index;

          V_rhs_value[bp_j] = 0.;
          V_rhs_index[bp_j] = index_j;

          for(int o_j=0; o_j<delta_j; o_j++){

            int j = o_j + bp_j*delta_j;

            overlap_indices& overlap_j = overlap_rhs[j];
            assert(index_j == overlap_j.index);

            //V_rhs_value[bp_j] += overlap_j.sign*psi_0(overlap_j.rhs, l0)*conj_value(psi_1(overlap_j.lhs, l1));
            V_rhs_value[bp_j] += overlap_j.sign*psi_0_vec[overlap_j.rhs]*conj_value(psi_1_vec[overlap_j.lhs]);
          }
        }
      }

      {
        for(size_t bp_j=0; bp_j<bp_rhs.size()-1; bp_j++){

          int index_j = V_rhs_index[bp_j];

          for(size_t bp_i=0; bp_i<bp_lhs.size()-1; bp_i++){

            int index_i = V_lhs_index[bp_i];

            overlap_r_r(index_i, index_j) += V_lhs_value[bp_i]*V_rhs_value[bp_j];
          }
        }
      }

      for(int r_j=0; r_j<r_dmn::dmn_size(); r_j++){
        for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++){

          int delta_r = rj_minus_ri(r_i, r_j);

          for(int s_j=0; s_j<s_dmn::dmn_size(); s_j++)
            for(int b_j=0; b_j<b_dmn::dmn_size(); b_j++)
              for(int s_i=0; s_i<s_dmn::dmn_size(); s_i++)
                for(int b_i=0; b_i<b_dmn::dmn_size(); b_i++)
                  overlap(b_i, s_i, b_j, s_j, delta_r)
                    += overlap_r_r(b_i, s_i, r_i, b_j, s_j, r_j);
        }
      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_sp_Greens_function<parameter_type, b_dmn, s_dmn, r_dmn>::set_up_psi_domain(int n_0, int Sz_0,
                                                                                              int n_1, int Sz_1)
    {
      double beta = parameters.get_beta();

      int N_0 = n_occupation_states(n_0, Sz_0);
      int N_1 = n_occupation_states(n_1, Sz_1);

      psi_lhs_domain::get_elements().resize(0);
      psi_rhs_domain::get_elements().resize(0);

      for(int l0=0; l0<N_0; l0++){
        double E_n0 = eigen_energies(n_0, Sz_0)[l0];
        double w_e   = std::exp(-beta*E_n0);

        if(w_e > CUT_OFF)
          psi_lhs_domain::get_elements().push_back(l0);
      }

      for(int l1=0; l1<N_1; l1++)
        psi_rhs_domain::get_elements().push_back(l1);

      psi_lhs_domain::get_size() = psi_lhs_domain::get_elements().size();
      psi_rhs_domain::get_size() = psi_rhs_domain::get_elements().size();
    }
  }

}

#endif
