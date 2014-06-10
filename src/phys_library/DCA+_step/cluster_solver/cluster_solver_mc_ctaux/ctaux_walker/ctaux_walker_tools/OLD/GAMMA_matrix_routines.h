//-*-C++-*-

#ifndef DCA_QMCI_GAMMA_TOOLS_H
#define DCA_QMCI_GAMMA_TOOLS_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \class   GAMMA_TOOLS
     *  \ingroup CT-AUX-WALKER
     *
     *  \author Peter Staar
     *  \brief  This class organizes the incremental construction of the \f$\Gamma\f$-matrix.
     *
     *  \f{eqnarray}{
     *   \Gamma_{i,j} &=& G_{i,j} - \delta_{i,j} * \frac{1+\gamma_j}{\gamma_j}                   (eqn 28) \\
     *   \gamma_{i,j} &=& e^{-\gamma(nu_j, mu_j) * HS_{field} * (HS_{spin}_new-HS_spin_old)} - 1 (eqn 17)
     *  \f}
     *
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    class GAMMA_TOOLS : public GAMMA_MATRIX_TOOLS<device_t>
    {
      typedef vertex_singleton        vertex_singleton_type;

      typedef typename parameters_type::Concurrency_Type concurrency_type;
      typedef typename parameters_type::profiler_type    profiler_t;
      
      typedef G_TOOLS<device_t, parameters_type> G_TOOLS_type;

    public:

      GAMMA_TOOLS(int                                                    id,
                  parameters_type&                                       parameters,
                  G_TOOLS<device_t, parameters_type>& G_TOOLS_obj_ref,
                  CV<parameters_type>&                CV_obj_ref);

      ~GAMMA_TOOLS();

      void initialize();

      template<class configuration_type>
      double compute_determinant_ratio_via_Gamma_LU(int                                configuration_e_spin_index,
                                                    HS_spin_states_type                new_HS_spin_value,
                                                    LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                    LIN_ALG::matrix<double, device_t>& N,
                                                    LIN_ALG::matrix<double, device_t>& G,
                                                    configuration_type&                configuration,
                                                    bool                               Bennett,
                                                    e_spin_states_type                 e_spin);

      template<class configuration_type>
      void compute_Gamma(LIN_ALG::matrix<double, device_t>& Gamma,
                         LIN_ALG::matrix<double, device_t>& N,
                         LIN_ALG::matrix<double, device_t>& G_precomputed,
                         configuration_type&                full_configuration,
                         e_spin_states_type                 e_spin);

      template<class configuration_type>
      void check_Gamma_LU(LIN_ALG::matrix<double, device_t>& Gamma_LU,
                          LIN_ALG::matrix<double, device_t>& N_up,
                          LIN_ALG::matrix<double, device_t>& G_precomputed,
                          configuration_type& configuration,
                          e_spin_states_type e_spin);

    private:

      using GAMMA_MATRIX_TOOLS<device_t>::initialize;
      using GAMMA_MATRIX_TOOLS<device_t>::resize;

      template<class configuration_type>
      double compute_Gamma_matrix_element(size_t                             vertex_index_1,
                                          size_t                             vertex_index_2,
                                          HS_spin_states_type                new_HS_spin_value,
                                          configuration_type&                configuration,
                                          LIN_ALG::matrix<double, device_t>& N,
                                          LIN_ALG::matrix<double, device_t>& G_precomputed,
                                          e_spin_states_type                 e_spin);

      template<class configuration_type>
      void apply_Bennett(int                                index,
                         configuration_type&                configuration,
                         LIN_ALG::matrix<double, device_t>& Gamma_LU,
                         LIN_ALG::matrix<double, device_t>& N,
                         LIN_ALG::matrix<double, device_t>& G_precomputed,
                         e_spin_states_type                 e_spin);

      template<class configuration_type>
      void compute_row_and_col_Gamma(int                                col_index,
                                     configuration_type&                full_configuration,
                                     LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                     LIN_ALG::matrix<double, device_t>& N,
                                     LIN_ALG::matrix<double, device_t>& G_precomputed,
                                     e_spin_states_type                 e_spin);

      template<class configuration_type>
      void compute_row_and_col_Gamma_simply(int                                col_index,
                                            configuration_type&                full_configuration,
                                            LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                            LIN_ALG::matrix<double, device_t>& N,
                                            LIN_ALG::matrix<double, device_t>& G_precomputed,
                                            e_spin_states_type                 e_spin);


      template<class configuration_type>
      void grow_Gamma_LU_iteratively(int                        configuration_e_spin_index,
                                     HS_spin_states_type        new_HS_spin_value,
                                     LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                     LIN_ALG::matrix<double, device_t>& N,
                                     LIN_ALG::matrix<double, device_t>& G_precomputed,
                                     configuration_type&        configuration,
                                     e_spin_states_type         e_spin);

      template<class configuration_type>
      void compute_s_and_w(int                                 j_index,
                           LIN_ALG::matrix<double, device_t>&  N,
                           LIN_ALG::matrix<double, device_t>&  Gamma_LU,
                           LIN_ALG::matrix<double, device_t>&  G_precomputed,
                           configuration_type&                 full_configuration,
                           e_spin_states_type                 e_spin);

      template<class configuration_type>
      void compute_s_and_w_2(int                                 j_index,
                             LIN_ALG::matrix<double, device_t>&  N,
                             LIN_ALG::matrix<double, device_t>&  Gamma_LU,
                             LIN_ALG::matrix<double, device_t>&  G_precomputed,
                             configuration_type&                 full_configuration,
                             e_spin_states_type                 e_spin);

      void solve_Gamma_L(LIN_ALG::matrix<double, device_t>& Gamma_LU);
      void solve_Gamma_U(LIN_ALG::matrix<double, device_t>& Gamma_LU);


      void print_Gamma_from_Gamma_LU(LIN_ALG::matrix<double, device_t>& Gamma_LU);

      template<class configuration_type>
      void print_Gamma(LIN_ALG::matrix<double, device_t>& Gamma_LU,
                       LIN_ALG::matrix<double, device_t>& N,
                       LIN_ALG::matrix<double, device_t>& G_precomputed,
                       configuration_type&                full_configuration,
                       e_spin_states_type                 e_spin);

    private:

      int thread_id;

      parameters_type&  parameters;
      concurrency_type& concurrency;

      CV<parameters_type>&                CV_obj;
      G_TOOLS<device_t, parameters_type>& G_TOOLS_obj;

      using GAMMA_MATRIX_TOOLS<device_t>::indices_up;
      using GAMMA_MATRIX_TOOLS<device_t>::indices_dn;

      using GAMMA_MATRIX_TOOLS<device_t>::exp_V_up;
      using GAMMA_MATRIX_TOOLS<device_t>::exp_V_dn;


      using GAMMA_MATRIX_TOOLS<device_t>::i_indices;
      using GAMMA_MATRIX_TOOLS<device_t>::j_indices;

      using GAMMA_MATRIX_TOOLS<device_t>::is_Bennett;

      using GAMMA_MATRIX_TOOLS<device_t>::exp_Vi;
      using GAMMA_MATRIX_TOOLS<device_t>::exp_Vj;

      using GAMMA_MATRIX_TOOLS<device_t>::r;
      using GAMMA_MATRIX_TOOLS<device_t>::c;
      using GAMMA_MATRIX_TOOLS<device_t>::delta;
    };

    template<LIN_ALG::device_type device_t, typename parameters_type>
    GAMMA_TOOLS<device_t, parameters_type>::GAMMA_TOOLS(int                                                    id,
                                                        parameters_type&                                       parameters_ref,
                                                        G_TOOLS<device_t, parameters_type>& G_TOOLS_obj_ref,
                                                        CV<parameters_type>&                CV_obj_ref):
      GAMMA_MATRIX_TOOLS<device_t>(id, parameters_ref.get_K_PHANI()),

      thread_id(id),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CV_obj     (CV_obj_ref),
      G_TOOLS_obj(G_TOOLS_obj_ref)
    {}

    template<LIN_ALG::device_type device_t, typename parameters_type>
    GAMMA_TOOLS<device_t, parameters_type>::~GAMMA_TOOLS()
    {}

    template<LIN_ALG::device_type device_t, typename parameters_type>
    void GAMMA_TOOLS<device_t, parameters_type>::initialize()
    {
      GAMMA_MATRIX_TOOLS<device_t>::initialize();
    }

    /*!
     *  \f{eqnarray}{
     *   \Gamma_{i,j} = G_{i,j} - \delta_{i,j} * (1+\gamma_j)/\gamma_j                    (eqn 28) \\
     *   \gamma_ij = \exp^{-gamma(nu_j, mu_j) * HS_field * (HS_spin_new-HS_spin_old)} - 1 (eqn 17)
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    double GAMMA_TOOLS<device_t, parameters_type>::compute_Gamma_matrix_element(size_t                             configuration_e_spin_index_i,
                                                                                size_t                             configuration_e_spin_index_j,
                                                                                HS_spin_states_type                new_HS_spin,
                                                                                configuration_type&                full_configuration,
                                                                                LIN_ALG::matrix<double, device_t>& N,
                                                                                LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                                e_spin_states_type                 e_spin)
    {
      std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

      assert(configuration_e_spin_index_i < size_t(configuration_e_spin.size() ));
      assert(configuration_e_spin_index_j < size_t(configuration_e_spin.size() ));

      vertex_singleton_type& v_i = configuration_e_spin[configuration_e_spin_index_i];
      vertex_singleton_type& v_j = configuration_e_spin[configuration_e_spin_index_j];

      int vertex_index_i = v_i.get_configuration_index();
      int vertex_index_j = v_j.get_configuration_index();

      double result;

      if(full_configuration[vertex_index_i].is_Bennett() || full_configuration[vertex_index_j].is_Bennett()){
        result = configuration_e_spin_index_i == configuration_e_spin_index_j ? 1. : 0. ;
      }
      else{
        result = G_TOOLS_obj.compute_G_matrix_element(configuration_e_spin_index_i, configuration_e_spin_index_j, N, G_precomputed, configuration_e_spin);

        if(configuration_e_spin_index_i == configuration_e_spin_index_j){
          double gamma_k = CV_obj.exp_delta_V(v_j, new_HS_spin);
          result -= (gamma_k)/(gamma_k-1.);
        }
      }

      return result;
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::compute_Gamma(LIN_ALG::matrix<double, device_t>& Gamma,
                                                               LIN_ALG::matrix<double, device_t>& N,
                                                               LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                               configuration_type&                full_configuration,
                                                               e_spin_states_type                 e_spin)
    {
      Gamma.resize(full_configuration.get_changed_spin_indices_e_spin(e_spin).size());

      assert(Gamma.get_number_of_rows() == Gamma.get_number_of_cols());

      for(int i=0; i<Gamma.get_number_of_rows(); i++){
        for(int j=0; j<Gamma.get_number_of_cols(); j++){

          int v_ind_i = full_configuration.get_changed_spin_indices_e_spin(e_spin)[i];
          int v_ind_j = full_configuration.get_changed_spin_indices_e_spin(e_spin)[j];

          HS_spin_states_type new_HS_i = full_configuration.get_changed_spin_values_e_spin(e_spin)[i];

          Gamma(i,j) = compute_Gamma_matrix_element(v_ind_i, v_ind_j, new_HS_i, full_configuration, N, G_precomputed, e_spin);
        }
      }
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    double GAMMA_TOOLS<device_t, parameters_type>::compute_determinant_ratio_via_Gamma_LU(int                                configuration_e_spin_index,
                                                                                          HS_spin_states_type                new_HS_spin,
                                                                                          LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                                                          LIN_ALG::matrix<double, device_t>& N,
                                                                                          LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                                          configuration_type&                full_configuration,
                                                                                          bool                               Bennett,
                                                                                          e_spin_states_type                 e_spin)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, "CT-AUX", __LINE__);

      double determinant_ratio;

      std::vector<int>&                   changed_spin_indices = full_configuration.get_changed_spin_indices_e_spin(e_spin);
      std::vector<HS_spin_states_type>&   changed_spin_values  = full_configuration.get_changed_spin_values_e_spin(e_spin);
      std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

      GAMMA_MATRIX_TOOLS<device_t>::resize(changed_spin_indices.size(), e_spin);

      if(Bennett)
        {
          //profiler_t profiler(concurrency, "Bennet-case", __FUNCTION__, __LINE__);

          int changed_index = find(changed_spin_indices.begin(), changed_spin_indices.end(), configuration_e_spin_index) - changed_spin_indices.begin();

          assert(changed_index < Gamma_LU.get_current_size().first && changed_index >= 0);

          determinant_ratio = 1./LIN_ALG::LU_MATRIX_OPERATIONS<device_t>::det(Gamma_LU);

          apply_Bennett(changed_index, full_configuration, Gamma_LU, N, G_precomputed, e_spin);

          determinant_ratio *= LIN_ALG::LU_MATRIX_OPERATIONS<device_t>::det(Gamma_LU);

          int spin_orbital                = configuration_e_spin[configuration_e_spin_index].get_spin_orbital();
          int spin_orbital_paired         = configuration_e_spin[configuration_e_spin_index].get_paired_spin_orbital();
          HS_spin_states_type old_HS_spin = changed_spin_values[changed_index];
          HS_field_sign_type  HS_field    = configuration_e_spin[configuration_e_spin_index].get_HS_field();

          assert(old_HS_spin != HS_ZERO && new_HS_spin == HS_ZERO);

          double phani_gamma = CV_obj.exp_delta_V(spin_orbital, spin_orbital_paired, old_HS_spin, new_HS_spin, HS_field)-1.;

          determinant_ratio = -determinant_ratio/phani_gamma;
        }
      else
        {
          //profiler_t profiler(concurrency, "general-case", __FUNCTION__, __LINE__);

          grow_Gamma_LU_iteratively(configuration_e_spin_index, new_HS_spin, Gamma_LU, N, G_precomputed, full_configuration, e_spin);

          vertex_singleton_type& v_new = configuration_e_spin[configuration_e_spin_index];
          double phani_gamma           = CV_obj.exp_delta_V(v_new, new_HS_spin)-1.;

          int n             = Gamma_LU.get_current_size().first;
          determinant_ratio = -phani_gamma*LIN_ALG::MEMORY_MANAGEMENT<device_t>::get(Gamma_LU.get_ptr(n-1, n-1));//Gamma_LU.get(n-1,n-1);
        }

      {
        //profiler_t profiler(concurrency, "LU-ratio", "Gamma-tools", __LINE__);

        double ratio = LIN_ALG::LU_MATRIX_OPERATIONS<device_t>::ratio(Gamma_LU);
        if(ratio > 1.e6)
          return 0;
      }

      return determinant_ratio;
    }

    /*!                 /            |   \          /            |      \ /            |      \
     *                  |            |   |          |            |      | |            |      |
     *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
     *                  |            |   |  ---/    |            |      | |            |      |
     *                  |------------|---|          |------------|---   | |------------|------|
     *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::grow_Gamma_LU_iteratively(int                                configuration_e_spin_index,
                                                                           HS_spin_states_type                new_HS_spin_value,
                                                                           LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                                           LIN_ALG::matrix<double, device_t>& N,
                                                                           LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                           configuration_type&                full_configuration,
                                                                           e_spin_states_type                 e_spin)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, "CT-AUX", __LINE__);

      Gamma_LU.add_row_and_col();
      int n = Gamma_LU.get_number_of_rows();

      compute_s_and_w(configuration_e_spin_index, N, Gamma_LU, G_precomputed, full_configuration, e_spin);

      solve_Gamma_L(Gamma_LU);
      solve_Gamma_U(Gamma_LU);

      {// compute \beta
        //profiler_t profiler(concurrency, __FUNCTION__, "Gamma-tools", __LINE__);

        vertex_singleton_type& v_j = full_configuration.get(e_spin)[configuration_e_spin_index];
        double gamma_k             = CV_obj.exp_delta_V(v_j, new_HS_spin_value);

        {
          double xy = LIN_ALG::DOT<device_t>::execute(n-1, Gamma_LU.get_ptr(n-1,0), Gamma_LU.get_leading_dimension(), Gamma_LU.get_ptr(0,n-1), 1);
          LIN_ALG::MEMORY_MANAGEMENT<device_t>::add(Gamma_LU.get_ptr(n-1, n-1), -(gamma_k)/(gamma_k-1.)-xy);
        }
      }
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::compute_s_and_w(int                                j_index,
                                                                 LIN_ALG::matrix<double, device_t>& N,
                                                                 LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                                 LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                 configuration_type&                full_configuration,
                                                                 e_spin_states_type                 e_spin)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, "Gamma-tools", __LINE__);

      int n = Gamma_LU.get_number_of_rows();

      std::vector<int>&                   changed_spin_indices = full_configuration.get_changed_spin_indices_e_spin(e_spin);
      std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

      LIN_ALG::vector<int   , device_t>& indices = GAMMA_MATRIX_TOOLS<device_t>::get_indices(e_spin);
      LIN_ALG::vector<double, device_t>& exp_V   = GAMMA_MATRIX_TOOLS<device_t>::get_exp_V  (e_spin);

      assert(indices.difference(full_configuration.get_changed_spin_indices_e_spin(e_spin))<1.e-6);

      GAMMA_MATRIX_TOOLS<device_t>::push_back_new_values(j_index, CV_obj.exp_V(configuration_e_spin[j_index]), e_spin);

      G_TOOLS_obj.compute_row_on_Gamma_matrix(indices.size()-1, indices, exp_V, N, G_precomputed, Gamma_LU.get_ptr(n-1, 0  ), Gamma_LU.get_leading_dimension());
      G_TOOLS_obj.compute_col_on_Gamma_matrix(indices.size()-1, indices, exp_V, N, G_precomputed, Gamma_LU.get_ptr(0  , n-1), 1.);

      for(int l=0; l<int(changed_spin_indices.size()); ++l){

        int  vertex_index_i = configuration_e_spin[changed_spin_indices[l]].get_configuration_index();
        bool is_Bennett     = full_configuration[vertex_index_i].is_Bennett();

        if(is_Bennett){
          LIN_ALG::MEMORY_MANAGEMENT<device_t>::set(Gamma_LU.get_ptr(n-1, l  ), 0.);
          LIN_ALG::MEMORY_MANAGEMENT<device_t>::set(Gamma_LU.get_ptr(l  , n-1), 0.);

          LIN_ALG::MEMORY_MANAGEMENT<device_t>::set(Gamma_LU.get_ptr(l,l), 1.);
        }
      }
    }

    /*!
     *  \f{eqnarray}{
     *   y &=& L^{-1} s
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    void GAMMA_TOOLS<device_t, parameters_type>::solve_Gamma_L(LIN_ALG::matrix<double, device_t>& Gamma_LU)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, "Gamma-tools", __LINE__);

      int n   = Gamma_LU.get_current_size().first;
      int lda = Gamma_LU.get_global_size().first;

      LIN_ALG::TRSV<device_t>::execute('L', 'N', 'U', n-1, Gamma_LU.get_ptr(0,0), lda, Gamma_LU.get_ptr(0,n-1), 1);
    }

    /*!
     *  \f{eqnarray}{
     *    x^T &=& U^{-1} w
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    void GAMMA_TOOLS<device_t, parameters_type>::solve_Gamma_U(LIN_ALG::matrix<double, device_t>& Gamma_LU)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, "Gamma-tools", __LINE__);

      int n   = Gamma_LU.get_current_size().first;
      int lda = Gamma_LU.get_global_size().first;

      LIN_ALG::TRSV<device_t>::execute('U', 'T', 'N', n-1, Gamma_LU.get_ptr(0,0), lda, Gamma_LU.get_ptr(n-1,0), lda);
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::apply_Bennett(int                                changed_index,
                                                               configuration_type&                full_configuration,
                                                               LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                               LIN_ALG::matrix<double, device_t>& N,
                                                               LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                               e_spin_states_type                 e_spin)
    {
      int n  = Gamma_LU.get_number_of_rows();
      int LD = Gamma_LU.get_leading_dimension();

      r    .resize(n);
      c    .resize(n);
      delta.resize(n);

      compute_row_and_col_Gamma(changed_index, full_configuration, Gamma_LU, N, G_precomputed, e_spin);

      {
        GAMMA_MATRIX_TOOLS<device_t>::set_delta(changed_index, -1.);
        LIN_ALG::BENNET<device_t>::standard_Bennet(n, LD, Gamma_LU.get_ptr(), c.get_ptr(), delta.get_ptr());
      }

      {
        GAMMA_MATRIX_TOOLS<device_t>::set_delta(changed_index, -1.);
        LIN_ALG::BENNET<device_t>::standard_Bennet(n, LD, Gamma_LU.get_ptr(), delta.get_ptr(), r.get_ptr());
      }
    }


    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::compute_row_and_col_Gamma(int                                col_index,
                                                                           configuration_type&                full_configuration,
                                                                           LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                                           LIN_ALG::matrix<double, device_t>& N,
                                                                           LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                           e_spin_states_type                 e_spin)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, "Gamma-tools", __LINE__);

      std::vector<int>&                   changed_spin_indices = full_configuration.get_changed_spin_indices_e_spin(e_spin);
      std::vector<HS_spin_states_type>&   changed_spin_values  = full_configuration.get_changed_spin_values_e_spin(e_spin);
      std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

      int changed_index_e_spin = changed_spin_indices[col_index];

      LIN_ALG::vector<int   , device_t>& indices = GAMMA_MATRIX_TOOLS<device_t>::get_indices(e_spin);
      LIN_ALG::vector<double, device_t>& exp_V   = GAMMA_MATRIX_TOOLS<device_t>::get_exp_V  (e_spin);

      assert(indices.difference(full_configuration.get_changed_spin_indices_e_spin(e_spin))<1.e-6);

      G_TOOLS_obj.compute_col_on_Gamma_matrix(col_index, indices, exp_V, N, G_precomputed, c.get_ptr(), 1);
      G_TOOLS_obj.compute_row_on_Gamma_matrix(col_index, indices, exp_V, N, G_precomputed, r.get_ptr(), 1);

      for(int l=0; l<int(changed_spin_indices.size()); ++l){
        //int  vertex_index_i = configuration_e_spin[indices[l]].get_configuration_index();
        int  vertex_index_i = configuration_e_spin[changed_spin_indices[l]].get_configuration_index();
        bool is_Bennett     = full_configuration[vertex_index_i].is_Bennett();

        if(is_Bennett){
          LIN_ALG::MEMORY_MANAGEMENT<device_t>::set(r.get_ptr(l), 0.);
          LIN_ALG::MEMORY_MANAGEMENT<device_t>::set(c.get_ptr(l), 0.);
        }
      }

      {
        HS_spin_states_type current_HS_spin = changed_spin_values[col_index];
        double gamma_k                      = CV_obj.exp_delta_V(configuration_e_spin[changed_index_e_spin], current_HS_spin);

        LIN_ALG::MEMORY_MANAGEMENT<device_t>::set(r.get_ptr(col_index), 0.);
        LIN_ALG::MEMORY_MANAGEMENT<device_t>::add(c.get_ptr(col_index), -(gamma_k)/(gamma_k-1.)-1.); // notice that the extra minus 1 is there to set diagonal to unity !
      }
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::compute_row_and_col_Gamma_simply(int                                col_index,
                                                                                  configuration_type&                full_configuration,
                                                                                  LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                                                  LIN_ALG::matrix<double, device_t>& N,
                                                                                  LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                                  e_spin_states_type                 e_spin)

    {
      int n = Gamma_LU.get_number_of_rows();

      std::vector<int>&                 changed_spin_indices = full_configuration.get_changed_spin_indices_e_spin(e_spin);
      std::vector<HS_spin_states_type>& changed_spin_values  = full_configuration.get_changed_spin_values_e_spin (e_spin);

      int changed_index_e_spin = changed_spin_indices[col_index];

      for(int i=0; i<n; i++){

        int index_e_spin_i                  = changed_spin_indices[i];
        HS_spin_states_type current_HS_spin = changed_spin_values [i];

        c[i] = compute_Gamma_matrix_element(index_e_spin_i, changed_index_e_spin, current_HS_spin, full_configuration, N, G_precomputed, e_spin);
        r[i] = compute_Gamma_matrix_element(changed_index_e_spin, index_e_spin_i, current_HS_spin, full_configuration, N, G_precomputed, e_spin);
      }

      c[col_index] += -1.;
      r[col_index]  =  0.;
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::check_Gamma_LU(LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                                LIN_ALG::matrix<double, device_t>& N,
                                                                LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                configuration_type&                full_configuration,
                                                                e_spin_states_type                 e_spin)
    {
      assert(Gamma_LU.get_number_of_rows() == Gamma_LU.get_number_of_cols());

      double diff=0;

      for(int i=0; i<Gamma_LU.get_number_of_rows(); i++){
        for(int j=0; j<Gamma_LU.get_number_of_cols(); j++){

          int v_ind_i = full_configuration.get_changed_spin_indices_e_spin(e_spin)[i];
          int v_ind_j = full_configuration.get_changed_spin_indices_e_spin(e_spin)[j];

          HS_spin_states_type new_HS_i = full_configuration.get_changed_spin_values_e_spin(e_spin)[i];

          double Gamma_i_j = compute_Gamma_matrix_element(v_ind_i, v_ind_j, new_HS_i, full_configuration, N, G_precomputed, e_spin);

          double other  = 0;
          double factor = 1;

          for(int l=0; l<Gamma_LU.get_number_of_rows(); l++)
            {
              factor = 1;

              if(l > i || l > j)
                factor = 0;

              if(l == i && !(l > j) )
                factor = 1/Gamma_LU(i,i);

              other += Gamma_LU(i,l)*Gamma_LU(l,j) * factor;
            }

          if(fabs(Gamma_i_j-other) > diff ){
            diff = fabs(Gamma_i_j-other);
          }
        }
      }

      if(diff>1.e-6)
        {
          print_Gamma_from_Gamma_LU(Gamma_LU);
          //print_Gamma(Gamma_LU, N, G_precomputed, full_configuration, e_spin);

          for(int i=0; i<Gamma_LU.get_number_of_rows(); i++){
            cout << "\t" << i << "\t";
            for(int j=0; j<Gamma_LU.get_number_of_cols(); j++){

              int v_ind_i = full_configuration.get_changed_spin_indices_e_spin(e_spin)[i];
              int v_ind_j = full_configuration.get_changed_spin_indices_e_spin(e_spin)[j];

              HS_spin_states_type new_HS_i = full_configuration.get_changed_spin_values_e_spin(e_spin)[i];

              double Gamma_i_j = compute_Gamma_matrix_element(v_ind_i, v_ind_j, new_HS_i, full_configuration, N, G_precomputed, e_spin);

              double other  = 0;
              double factor = 1;

              for(int l=0; l<Gamma_LU.get_number_of_rows(); l++)
                {
                  factor = 1;

                  if(l > i || l > j)
                    factor = 0;

                  if(l == i && !(l > j) )
                    factor = 1/Gamma_LU(i,i);

                  other += Gamma_LU(i,l)*Gamma_LU(l,j) * factor;
                }

              if(fabs(Gamma_i_j-other)<1.e-6)
                cout << 0. << "\t";
              else
                cout << fabs(Gamma_i_j-other) << "\t";
            }
            cout << endl;
          }
          cout << endl;

          throw std::logic_error(__FUNCTION__);
        }
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void GAMMA_TOOLS<device_t, parameters_type>::print_Gamma(LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                                             LIN_ALG::matrix<double, device_t>& N,
                                                             LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                             configuration_type&                full_configuration,
                                                             e_spin_states_type                 e_spin)
    {
      cout << __FUNCTION__ << endl << endl;

      cout << scientific;
      cout.precision(6);

      for(int i=0; i<Gamma_LU.get_number_of_rows(); i++){
        for(int j=0; j<Gamma_LU.get_number_of_cols(); j++){

          int v_ind_i = full_configuration.get_changed_spin_indices_e_spin(e_spin)[i];
          int v_ind_j = full_configuration.get_changed_spin_indices_e_spin(e_spin)[j];

          HS_spin_states_type new_HS_i = full_configuration.get_changed_spin_values_e_spin(e_spin)[i];

          double Gamma_i_j = compute_Gamma_matrix_element(v_ind_i, v_ind_j, new_HS_i, full_configuration, N, G_precomputed, e_spin);
          cout << "\t" << Gamma_i_j;
        }
        cout << "\n";
      }
      cout << "\n";
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    void GAMMA_TOOLS<device_t, parameters_type>::print_Gamma_from_Gamma_LU(LIN_ALG::matrix<double, device_t>& Gamma_LU)
    {
      cout << __FUNCTION__ << endl << endl;

      cout << scientific;
      cout.precision(6);

      for(int i=0; i<Gamma_LU.get_number_of_rows(); i++){
        cout << "\t" << i << "\t";
        for(int j=0; j<Gamma_LU.get_number_of_cols(); j++){

          double other  = 0;
          double factor = 1;

          for(int l=0; l<Gamma_LU.get_number_of_rows(); l++)
            {
              factor = 1;

              if(l > i || l > j)
                factor = 0;

              if(l == i && !(l > j) )
                factor = 1/Gamma_LU(i,i);

              other += Gamma_LU(i,l)*Gamma_LU(l,j) * factor;
            }

          cout << "\t" << other;
        }
        cout << endl;
      }

      cout << endl;
    }

  }

}

#endif

