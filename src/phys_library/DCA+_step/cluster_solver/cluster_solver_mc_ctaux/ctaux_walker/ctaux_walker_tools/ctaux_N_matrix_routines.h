//-*-C++-*-

#ifndef DCA_QMCI_N_MATRIX_ROUTINES_H
#define DCA_QMCI_N_MATRIX_ROUTINES_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \class   N_TOOLS
     *  \ingroup CT-AUX-WALKER
     *
     *  \author Peter Staar
     *  \brief  This class organizes the construction of the N-matrix.
     *
     *  The N-matrix can be computed directly from the Hubbard-spin configuration,
     *
     *  \f{eqnarray}{
     *     N^{-1} = e^{-\gamma \: HS_{field} \: HS_{spin}} - G^{0} * [e^{- \gamma * HS_{field} * HS_{spin}}-1]
     *  \f}
     *
     *  If non-intracting spins are added to the configuration, then the N-matrix needs to be updated (--> see formula 46 ),
     *
     *  \f{eqnarray}{
     *    N_{i,j} &=& N_{i,j} \mbox{ if } i<n && j<n
     *    N_{i,j} &=& \sum_{k} G_0_{i,k}*(exp(-\gamma*e_spin*HS_spin(k))-1*)*N_{k,j} \mbox{ if } i \leq n && j<n
     *    N_{i,j} &=& \delta_{i,j}  \mbox{ if } j \leq n
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    class N_TOOLS : public N_MATRIX_TOOLS<device_t, parameters_type>
    {
      const static int MAX_VERTEX_SINGLETS=4;

      typedef vertex_singleton    vertex_singleton_type;

      typedef typename parameters_type::concurrency_type concurrency_type;
      typedef typename parameters_type::profiler_type    profiler_t;

    public:

      N_TOOLS(int                                     id,
              parameters_type&                        parameters,
              CV<parameters_type>& CV_obj_ref);

      ~N_TOOLS();

      double get_Gflop();

      template<class configuration_type>
      void build_N_matrix(configuration_type&                configuration,
                          LIN_ALG::matrix<double, device_t>& N,
                          LIN_ALG::matrix<double, device_t>& G0,
                          e_spin_states_type                 e_spin);

      template<class configuration_type>
      void update_N_matrix(configuration_type&                full_configuration,
                           LIN_ALG::matrix<double, device_t>& G0,
                           LIN_ALG::matrix<double, device_t>& N,
                           e_spin_states_type                 e_spin);

      template<class configuration_type>
      void rebuild_N_matrix_via_Gamma_LU(configuration_type&                full_configuration,
                                         LIN_ALG::matrix<double, device_t>& N,
                                         LIN_ALG::matrix<double, device_t>& Gamma_LU,
                                         LIN_ALG::matrix<double, device_t>& G,
                                         e_spin_states_type                 e_spin);

      template<class configuration_type>
      void check_N_matrix(configuration_type&        configuration,
                          LIN_ALG::matrix<double, device_t>& G0,
                          LIN_ALG::matrix<double, device_t>& N,
                          LIN_ALG::matrix<double, device_t>& Gamma,
                          e_spin_states_type         e_spin);

    private:

      void copy_rows_N(std::vector<int>&                  permutation,
                       LIN_ALG::matrix<double, device_t>& N);

      template<class configuration_type>
      void compute_G_changed_vertices_to_all_vertex(LIN_ALG::matrix<double, device_t>& N,
                                                    LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                    configuration_type&                full_configuration,
                                                    e_spin_states_type                 e_spin);

      void compute_d_vector(std::vector<int>&                   permutation,
                            LIN_ALG::matrix<double, device_t>&  N,
                            std::vector<HS_spin_states_type>&   spin_values,
                            std::vector<vertex_singleton_type>& configuration_e_spin,
                            LIN_ALG::vector<double, LIN_ALG::CPU>&                d_inv);

      void scale_rows_N(std::vector<int>&                  permutation,
                        LIN_ALG::vector<double, LIN_ALG::CPU>&  d_inv,
                        LIN_ALG::matrix<double, device_t>& N);

      template<class configuration_type>
      static bool assert_N_matrix_format(configuration_type&  full_configuration);

      template<class configuration_type>
      static bool assert_that_there_are_no_Bennett_spins(configuration_type& full_configuration);

    private:

      int thread_id;
      int stream_id;

      double GFLOP;

      parameters_type&  parameters;
      concurrency_type& concurrency;

      CV<parameters_type>& CV_obj;

      LIN_ALG::vector<double, LIN_ALG::CPU> exp_gamma_s, one_min_exp_gamma_s, d_inv, exp_V_minus_one_val;

      LIN_ALG::matrix<double, device_t> G;
      LIN_ALG::matrix<double, device_t> N_new_spins;
      LIN_ALG::matrix<double, device_t> G0_times_exp_V_minus_one;
    };

    template<LIN_ALG::device_type device_t, typename parameters_type>
    N_TOOLS<device_t, parameters_type>::N_TOOLS(int                                     id,
                                                parameters_type&                        parameters_ref,
                                                CV<parameters_type>& CV_obj_ref):

      N_MATRIX_TOOLS<device_t, parameters_type>(id, parameters_ref),

      thread_id(id),
      stream_id(0),

      GFLOP(0.),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CV_obj(CV_obj_ref),

      exp_gamma_s             ("exp_gamma_s (N_TOOLS)"        , MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      one_min_exp_gamma_s     ("one_min_exp_gamma_s (N_TOOLS)", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),

      d_inv                   ("d_inv (N_TOOLS)"              , MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      exp_V_minus_one_val     ("exp_V_minus_one_val (N_TOOLS)", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),

      //     G                       ("G (N_TOOLS)"                       , MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      //     N_new_spins             ("N_new_spins (N_TOOLS)"             , MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      //     G0_times_exp_V_minus_one("G0_times_exp_V_minus_one (N_TOOLS)", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI())

      G                       ("G (N_TOOLS)"                       , std::pair<int,int>(0,0), std::pair<int,int>(parameters.get_initial_matrix_size(), MAX_VERTEX_SINGLETS*parameters.get_K_PHANI())),
      N_new_spins             ("N_new_spins (N_TOOLS)"             , std::pair<int,int>(0,0), std::pair<int,int>(MAX_VERTEX_SINGLETS*parameters.get_K_PHANI(), parameters.get_initial_matrix_size())),
      G0_times_exp_V_minus_one("G0_times_exp_V_minus_one (N_TOOLS)", std::pair<int,int>(0,0), std::pair<int,int>(MAX_VERTEX_SINGLETS*parameters.get_K_PHANI(), parameters.get_initial_matrix_size()))
    {}

    template<LIN_ALG::device_type device_t, typename parameters_type>
    N_TOOLS<device_t, parameters_type>::~N_TOOLS()
    {}

    template<LIN_ALG::device_type device_t, typename parameters_type>
    double N_TOOLS<device_t, parameters_type>::get_Gflop()
    {
      double result = GFLOP;
      GFLOP         = 0.;

      if(result<0)
        std::cout << __FUNCTION__ << "\t" << result << "\n";

      return result;
    }

    /*!
     *  The N-matrix can be computed directly from the Hubbard-spin configuration,
     *
     *  \f{eqnarray}{
     *     N^{-1} = e^{-\gamma \: HS_{field} \: HS_{spin}} - G^{0} * [e^{- \gamma * HS_{field} * HS_{spin}}-1]
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void N_TOOLS<device_t, parameters_type>::build_N_matrix(configuration_type&                configuration,
                                                            LIN_ALG::matrix<double, device_t>& N,
                                                            LIN_ALG::matrix<double, device_t>& G0,
                                                            e_spin_states_type                 e_spin)
    {
      std::vector<vertex_singleton_type>&  configuration_e_spin = configuration.get(e_spin);
      int configuration_size(configuration_e_spin.size());

      // All interaction pairs are of the same spin type, which leads to a zero configuration size for one of the spin types.
      if (configuration_size == 0) {
        return;
      }

      exp_gamma_s        .resize(configuration_size);
      one_min_exp_gamma_s.resize(configuration_size);

      N.resize_no_copy(configuration_size);

      for(int i=0; i<configuration_size; ++i)
        exp_gamma_s[i] = CV_obj.exp_V(configuration_e_spin[i]);

      { // GEMD
        for(int i=0; i<configuration_size; ++i)
          one_min_exp_gamma_s[i] = (1.-exp_gamma_s[i]);

        LIN_ALG::GEMD<device_t>::execute(G0, N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(one_min_exp_gamma_s), N, thread_id, stream_id);
      }

      {
        double* exp_gamma_s_ptr = N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(exp_gamma_s);

        LIN_ALG::AXPY<device_t>::execute(configuration_size, 1., exp_gamma_s_ptr,
                                         1, N.get_ptr(), N.get_leading_dimension()+1,
                                         thread_id, stream_id);
      }

      LIN_ALG::GEINV<device_t>::execute(N);
    }


    /*!
     *  If non-intracting spins are added to the configuration, then the N-matrix needs to be updated (--> see formula 46 ),
     *
     *  \f{eqnarray}{
     *    N_{i,j} &=& N_{i,j} \mbox{ if } i<n && j<n \\
     *    N_{i,j} &=& \sum_{k} G_0_{i,k}*(exp(-\gamma*e_spin*HS_spin(k))-1*)*N_{k,j} \mbox{ if } i \leq n && j<n \\
     *    N_{i,j} &=& \delta_{i,j}  \mbox{ if } j \leq n
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void  N_TOOLS<device_t, parameters_type>::update_N_matrix(configuration_type&                configuration,
                                                              LIN_ALG::matrix<double, device_t>& G0,
                                                              LIN_ALG::matrix<double, device_t>& N,
                                                              e_spin_states_type                 e_spin)
    {
      //profiler_t profiler(concurrency, "update_N_matrix", "CT-AUX", __LINE__, true);

      std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
      int                                 configuration_size   = configuration_e_spin.size();

      // All interaction pairs are of the same spin type, which leads to a zero configuration size for one of the spin types.
      if (configuration_size == 0) {
        return;
      }

      int first_non_interacting_vertex_index = configuration.get_first_non_interacting_spin_index(e_spin);
      int first_shuffled_vertex_index        = configuration.get_first_shuffled_spin_index       (e_spin);

      assert(configuration.get_changed_spin_indices().size()            == 0);
      assert(configuration.get_changed_spin_indices_e_spin(e_UP).size() == 0);
      assert(configuration.get_changed_spin_indices_e_spin(e_DN).size() == 0);
      assert(configuration.assert_block_form(e_spin));
      assert(first_shuffled_vertex_index >= first_non_interacting_vertex_index);

      if (first_non_interacting_vertex_index == configuration_size) {
        assert(configuration_size == N.get_current_size().first);
        // std::cout << __FUNCTION__
        //           << "\tconfiguration_size = " << configuration_size
        //           << "\tN.get_current_size().first = " << N.get_current_size().first << std::endl;
        return;
      }

      {
        N.resize(configuration_size);  // Move after the next if block?
      }
      
      { // set columns to unity ...
        //profiler_t profiler(concurrency, "(a) set columns to unity", __FUNCTION__, __LINE__, true);

        int i = first_non_interacting_vertex_index;

        int N_r = N.get_number_of_rows();
        int N_c = N.get_number_of_cols();

        int LD = N.get_leading_dimension();

        assert(N_r==N_c);
        LIN_ALG::LASET<device_t>::set_zero (    i, N_c-i, N.get_ptr(0,i), LD, thread_id, stream_id);
        LIN_ALG::LASET<device_t>::set_unity(N_r-i, N_c-i, N.get_ptr(i,i), LD, thread_id, stream_id);
      }

      if(first_shuffled_vertex_index == configuration_size || first_non_interacting_vertex_index == 0)
        return;

      {// G0_times_exp_V_minus_one <--- \sum_{l=0}^{l<vertex_index} G0_{i,l}*exp_V_minus_one[l]
        //profiler_t profiler(concurrency, "(b) GEMD", __FUNCTION__, __LINE__, true);

        std::pair<int,int> current_size(configuration_size-first_shuffled_vertex_index, first_non_interacting_vertex_index);

        G0_times_exp_V_minus_one.resize_no_copy(current_size);

        exp_V_minus_one_val.resize(first_non_interacting_vertex_index);
        for(int j=0; j<first_non_interacting_vertex_index; ++j)
          exp_V_minus_one_val[j] = CV_obj.exp_V(configuration_e_spin[j])-1.;

        double* diagonal_matrix_ptr = N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(exp_V_minus_one_val);

        LIN_ALG::GEMD<device_t>::execute(current_size,
                                         &G0.get_ptr()[first_shuffled_vertex_index], G0.get_leading_dimension(),
                                         diagonal_matrix_ptr,
                                         //N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(exp_V_minus_one_val),
                                         G0_times_exp_V_minus_one.get_ptr(),  G0_times_exp_V_minus_one.get_leading_dimension(),
                                         thread_id, stream_id);
      }

      { // G0_exp_V_minus_one * N
        //profiler_t profiler(concurrency, "(c) GEMM", __FUNCTION__, __LINE__, true);

        int m = configuration_size - first_shuffled_vertex_index;
        int k = first_non_interacting_vertex_index;
        int n = first_non_interacting_vertex_index;

        int LD_G0 = G0_times_exp_V_minus_one.get_global_size().first;
        int LD_N  = N.get_global_size().first;

        LIN_ALG::GEMM<device_t>::execute('N', 'N', m, n, k,
                                         1.,
                                         G0_times_exp_V_minus_one.get_ptr()       , LD_G0,
                                         N.get_ptr()                              , LD_N,
                                         0.,
                                         &N.get_ptr()[first_shuffled_vertex_index], LD_N,
                                         thread_id, stream_id);

        GFLOP += 2.*double(m)*double(k)*double(n)*(1.e-9);
      }
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void  N_TOOLS<device_t, parameters_type>::rebuild_N_matrix_via_Gamma_LU(configuration_type&                full_configuration,
                                                                            LIN_ALG::matrix<double, device_t>& N,
                                                                            LIN_ALG::matrix<double, device_t>& Gamma,
                                                                            LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                            e_spin_states_type                 e_spin)
    {
      //profiler_t profiler(concurrency, "rebuild_N_matrix_via_Gamma_LU", "CT-AUX", __LINE__, true);

      int Gamma_size = Gamma.get_current_size().first;

      if(Gamma_size == 0)
        return;

      std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
      int                                 configuration_size   = configuration_e_spin.size();  // What happens if configuration_size = 0?

      std::vector<HS_spin_states_type>& spin_values = full_configuration.get_changed_spin_values_e_spin(e_spin);
      std::vector<int>&                 permutation = full_configuration.get_changed_spin_indices_e_spin(e_spin);

      assert(spin_values.size()         ==     permutation.size());
      assert(Gamma_size                 == int(permutation.size()));
      assert(N.get_current_size().first == int(configuration_size));
      assert(assert_that_there_are_no_Bennett_spins(full_configuration));

      N_MATRIX_TOOLS<device_t, parameters_type>::set_permutation(permutation);

      { // get the rows of N corresponding to the new spins => N_new_spins
        //profiler_t profiler(concurrency, "(a) resize N && copy rows", __FUNCTION__, __LINE__, true);

        N_new_spins.resize_no_copy(std::pair<int,int>(Gamma_size, configuration_size));

        //copy_rows_N(permutation, N);
        N_MATRIX_TOOLS<device_t, parameters_type>::copy_rows(N, N_new_spins);
      }

      { // get the columns of G corresponding to the new spins => G_new_spins
        //profiler_t profiler(concurrency, "(b) resize G && copy cols", __FUNCTION__, __LINE__, true);

        G.resize_no_copy(std::pair<int,int>(configuration_size, Gamma_size));

        //if(true)
        {
          std::vector<double> exp_V(permutation.size());
          for(size_t l=0; l<permutation.size(); ++l)
            exp_V[l] = CV_obj.exp_V(configuration_e_spin[permutation[l]]);

          N_MATRIX_TOOLS<device_t, parameters_type>::compute_G_cols(exp_V, N, G_precomputed, G);
        }
        //else
        //compute_G_changed_vertices_to_all_vertex(N, G_precomputed, full_configuration, e_spin);
      }

      { // Gamma_LU * X = N(p_k,:) --> X = Gamma_inv_times_N_new_spins ==> (stored in N_new_spins)
        //profiler_t profiler(concurrency, "(c) LU-solve", __FUNCTION__, __LINE__, true);

        LIN_ALG::TRSM<device_t>::execute('L', 'U', Gamma, N_new_spins, thread_id, stream_id);
        LIN_ALG::TRSM<device_t>::execute('U', 'N', Gamma, N_new_spins, thread_id, stream_id);

        GFLOP += 2.*double(Gamma_size)*double(Gamma_size)*double(configuration_size)*(1.e-9);
      }

      { // do N - G*Gamma_inv_times_N_new_spins --> N  || DGEMM --> work-horsegg
        //profiler_t profiler(concurrency, "(d) dgemm", __FUNCTION__, __LINE__, true);

        LIN_ALG::GEMM<device_t>::execute(-1., G, N_new_spins, 1., N, thread_id, stream_id);

        GFLOP += 2.*double(configuration_size)*double(Gamma_size)*double(configuration_size)*(1.e-9);
      }

      { // do N*D_i --> N ( = final N !)
        //profiler_t profiler(concurrency, "(e) rescale", __FUNCTION__, __LINE__, true);

        compute_d_vector(permutation, N, spin_values, configuration_e_spin, d_inv);

        N_MATRIX_TOOLS<device_t, parameters_type>::scale_rows(N);
      }
    }

    /*
      template<LIN_ALG::device_type device_t, typename parameters_type>
      inline void N_TOOLS<device_t, parameters_type>::set_data()
      {
      std::vector<double> exp_V(permutation.size());
      std::vector<double> d_vec(permutation.size());

      {
      for(size_t l=0; l<permutation.size(); ++l)
      exp_V[l] = CV_obj.exp_V(configuration_e_spin[permutation[l]]);
      }

      {
      int                 spin_orbital, spin_orbital_paired;
      double              exp_delta_V;

      HS_field_sign       HS_field_sign;
      HS_spin_states_type old_HS_spin, new_HS_spin;

      for(int i=0; i<int(permutation.size()); ++i){

      HS_spin_states_type old_HS_spin = configuration_e_spin[permutation[i]].get_HS_spin();
      HS_spin_states_type new_HS_spin = spin_values[i];

      HS_field_sign HS_field_sign = configuration_e_spin[permutation[i]].get_HS_field();

      int spin_orbital        = configuration_e_spin[permutation[i]].get_spin_orbital();
      int spin_orbital_paired = configuration_e_spin[permutation[i]].get_paired_spin_orbital();

      if(old_HS_spin == HS_ZERO)
      {
      d_vec[i] = 1./CV_obj.exp_delta_V(spin_orbital, spin_orbital_paired, new_HS_spin, old_HS_spin, HS_field_sign);
      }
      else
      {
      d_inv[i] = 0;
      }
      }
      }

      N_MATRIX_TOOLS<device_t, parameters_type>::set_data(permutation, exp_V, d_vec);
      }
    */

    template<LIN_ALG::device_type device_t, typename parameters_type>
    inline void N_TOOLS<device_t, parameters_type>::copy_rows_N(std::vector<int>&                  permutation,
                                                                LIN_ALG::matrix<double, device_t>& N)
    {
      //profiler_t profiler(concurrency, __FUNCTION__, __FUNCTION__, __LINE__);

      if(true)
        {
          for(size_t i=0; i<permutation.size(); ++i)
            LIN_ALG::COPY<device_t>::row(N, permutation[i], N_new_spins, i);
        }
    }

    /*!
     *  If the vertex v_j was interacting, then
     *  \f{eqnarray}{
     *    G_{i,j} = (N_{i,j} e^{V_j} - delta_{i,j})/(e^{V_j} -1)
     *  \f}
     *
     *  If the vertex v_j was non-interacting, then
     *  \f{eqnarray}{
     *    G_{i,j} = \sum_l N_{i,l}*G0_{l,j}
     *  \f}
     */
    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void N_TOOLS<device_t, parameters_type>::compute_G_changed_vertices_to_all_vertex(LIN_ALG::matrix<double, device_t>& N,
                                                                                      LIN_ALG::matrix<double, device_t>& G_precomputed,
                                                                                      configuration_type&                full_configuration,
                                                                                      e_spin_states_type                 e_spin)
    {
      std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

      std::vector<int>&                   permutation = full_configuration.get_changed_spin_indices_e_spin(e_spin);

      int configuration_size = configuration_e_spin.size();
      int Gamma_size         = permutation.size();

      int N_ind = N.get_number_of_cols()-G_precomputed.get_number_of_cols();

      for(int l=0; l<Gamma_size; ++l){

        int j_ind = permutation[l];

        if(j_ind >= N_ind)
          LIN_ALG::COPY<device_t>::execute(configuration_size, G_precomputed.get_ptr(0,j_ind-N_ind), 1, G.get_ptr(0,l), 1);
        else
          {
            double exp_V       = CV_obj.exp_V(configuration_e_spin[j_ind]);
            double denumerator = exp_V-1.;

            double alpha = exp_V/denumerator;

            LIN_ALG::COPY <device_t>::execute(configuration_size, N.get_ptr(0,j_ind), 1, G.get_ptr(0,l), 1);
            LIN_ALG::SCALE<device_t>::execute(configuration_size, alpha, G.get_ptr(0,l), 1);

            //G(j_ind,l) -= 1./denumerator;
            LIN_ALG::MEMORY_MANAGEMENT<device_t>::add(G.get_ptr(j_ind,l), -1./denumerator);
          }
      }
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    inline void N_TOOLS<device_t, parameters_type>::compute_d_vector(std::vector<int>&                      permutation,
                                                                     LIN_ALG::matrix<double, device_t>&     /*N*/,
                                                                     std::vector<HS_spin_states_type>&      spin_values,
                                                                     std::vector<vertex_singleton_type>&    configuration_e_spin,
                                                                     LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv)
    {
      int                 spin_orbital, spin_orbital_paired;
      double              exp_delta_V;

      HS_field_sign       HS_field_sign;
      HS_spin_states_type old_HS_spin, new_HS_spin;

      d_inv.resize(permutation.size());

      std::vector<int> d_index(0);
      //std::vector<int> N_index(0);

      for(int i=0; i<int(permutation.size()); ++i){

        old_HS_spin   = configuration_e_spin[permutation[i]].get_HS_spin();
        new_HS_spin   = spin_values[i];

        HS_field_sign = configuration_e_spin[permutation[i]].get_HS_field();

        spin_orbital        = configuration_e_spin[permutation[i]].get_spin_orbital();
        spin_orbital_paired = configuration_e_spin[permutation[i]].get_paired_spin_orbital();

        if(old_HS_spin == HS_ZERO)
          {
            exp_delta_V = CV_obj.exp_delta_V(spin_orbital, spin_orbital_paired, new_HS_spin, old_HS_spin, HS_field_sign);
            d_inv[i] = 1./exp_delta_V;
          }
        else{
          d_inv[i] = 0;

          //d_inv[i] = 1./LIN_ALG::MEMORY_MANAGEMENT<device_t>::get(N.get_ptr(permutation[i],permutation[i]));
          //d_inv[i] = 1./N(permutation[i],permutation[i]);

          //d_index.push_back(i);
        }
      }

      //N_MATRIX_TOOLS<device_t, parameters_type>::set_d_vector(d_index, N, d_inv);

      N_MATRIX_TOOLS<device_t, parameters_type>::set_d_vector(d_inv);
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    inline void N_TOOLS<device_t, parameters_type>::scale_rows_N(std::vector<int>&                      permutation,
                                                                 LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv,
                                                                 LIN_ALG::matrix<double, device_t>&     N)
    {
      for(size_t i=0; i<permutation.size(); ++i)
        LIN_ALG::SCALE<device_t>::row(N, d_inv[i], permutation[i]);
    }












    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    bool N_TOOLS<device_t, parameters_type>::assert_that_there_are_no_Bennett_spins(configuration_type& full_configuration)
    {
      {
        std::vector<vertex_singleton_type>&    configuration_e_spin = full_configuration.get(e_UP);
        std::vector<int>&                 permutation = full_configuration.get_changed_spin_indices_e_spin(e_UP);

        for(size_t i=0; i< permutation.size() ; i++)
          {
            int configuration_index = configuration_e_spin[permutation[i]].get_configuration_index();
            //assert(!full_configuration[configuration_index].is_Bennett());
            if(full_configuration[configuration_index].is_Bennett()) throw std::logic_error(__FUNCTION__);
          }
      }

      {
        std::vector<vertex_singleton_type>&    configuration_e_spin = full_configuration.get(e_DN);
        std::vector<int>&                 permutation = full_configuration.get_changed_spin_indices_e_spin(e_DN);

        for(size_t i=0; i< permutation.size() ; i++)
          {
            int configuration_index = configuration_e_spin[permutation[i]].get_configuration_index();
            //assert(!full_configuration[configuration_index].is_Bennett());
            if(full_configuration[configuration_index].is_Bennett()) throw std::logic_error(__FUNCTION__);
          }
      }

      return true;
    }

    template<LIN_ALG::device_type device_t, typename parameters_type>
    template<class configuration_type>
    void N_TOOLS<device_t, parameters_type>::check_N_matrix(configuration_type&                configuration,
                                                            LIN_ALG::matrix<double, device_t>& N,
                                                            LIN_ALG::matrix<double, device_t>& G0,
                                                            LIN_ALG::matrix<double, device_t>& /*Gamma*/,
                                                            e_spin_states_type                 e_spin)
    {
      LIN_ALG::matrix<double, device_t> N_correct(N.get_current_size(), N.get_global_size());

      std::cout.precision(4);

      build_N_matrix(configuration, N_correct, G0, e_spin);

      N.difference(N_correct);

      //     if(! N.difference(N_correct))
      //       {
      //        configuration.print();
      //        configuration.print(e_spin);

      //        N_correct.print();
      //        N.print();

      //        std::cout << "\t\t e_spin : " << e_spin << std::endl;
      //        Gamma.print();

      //        for(int i=0; i<N.get_current_size(); i++)
      //          std::cout << "\t\t" << configuration[i].get_HS_spin();
      //        std::cout << "\n";

      //        //throw std::logic_error(__FUNCTION__);

      //       }
    }

  }

}

#endif








//   template<LIN_ALG::device_type device_t, typename parameters_type, class base_cluster_type>
//   template<class configuration_type, class vertex_vertex_matrix_type>
//   void N_TOOLS<device_t, parameters_type, base_cluster_type>::build_N_matrix(configuration_type&        configuration,
//                                                                 vertex_vertex_matrix_type& N,
//                                                                 vertex_vertex_matrix_type& G0,
//                                                                 e_spin_states_type         e_spin)
//   {
//     // N^-1 = exp(- gamma * HS_field * HS_spin) - G0 * [exp(- gamma * HS_field * HS_spin) -1]

//     std::vector<vertex_singleton_type>&  configuration_e_spin = configuration.get(e_spin);
//     int configuration_size(configuration_e_spin.size());

//     std::vector<double> exp_gamma_s(configuration_size);

//     N.resize_no_copy(configuration_size);

//     for(int i=0;i<configuration_size;++i){
//       HS_spin_states_type HS_spin  = configuration_e_spin[i].get_HS_spin();
//       HS_field_sign_type  HS_field = configuration_e_spin[i].get_HS_field();

//       int spin_orbital             = configuration_e_spin[i].get_spin_orbital();
//       int spin_orbital_paired      = configuration_e_spin[i].get_paired_spin_orbital();

//       double exp_V = CV_obj.exp_V(spin_orbital, spin_orbital_paired, HS_spin, HS_field);
//       exp_gamma_s[i]=exp_V;
//     }

//     { // GEMD
//       for(int i=0;i<configuration_size;++i)
//      for(int j=0;j<configuration_size;++j)
//        N(i,j) = -G0(i,j)*(exp_gamma_s[j]-1.);
//     }

//     { // axpy
//       for(int i=0;i<configuration_size;++i)
//      N(i,i)+=exp_gamma_s[i];
//     }

//     { // INVERT
//       invert_plan<double> invert_pln(N.get_current_size(), N.get_global_size());

//       memcpy(invert_pln.Matrix,  &(N(0,0)),         sizeof(double)*N.get_global_size()*N.get_global_size());
//       invert_pln.execute_plan();
//       memcpy(&(N(0,0)), invert_pln.inverted_matrix, sizeof(double)*N.get_global_size()*N.get_global_size());
//     }
//   }


//   template<LIN_ALG::device_type device_t, typename parameters_type, class base_cluster_type>
//   template<class configuration_type, class vertex_vertex_matrix_type>
//   void  N_TOOLS<device_t, parameters_type, base_cluster_type>::update_N_matrix(configuration_type&        configuration,
//                                                                             vertex_vertex_matrix_type& G0,
//                                                                             vertex_vertex_matrix_type& N,
//                                                                             e_spin_states_type         e_spin)
//   {
//     profiler_t profiler(concurrency, "TOOLS update_N_matrix", __FILE__, __LINE__);
//     // update the row of directly for all columns
//     // see formula (46)
//     // if configuration[i].get_HS_spin() == zero
//     // {
//     //    ifconfiguration[j].get_HS_spin() != zero
//     //      N_old->new --> N(i,j) --> \sum_{k \elem interacting} G_0(i,k)*(exp(-\gamma*e_spin*HS_spin(k))-1*)*N(k,j)
//     //    else
//     //      N_old->old --> N(i,j) --> \delta_{i,j}
//     // }
//     //

//     std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
//     int                                 configuration_size   = configuration_e_spin.size();

//     int first_non_interacting_vertex_index = configuration.get_first_non_interacting_spin_index(e_spin);
//     int first_shuffled_vertex_index        = configuration.get_first_shuffled_spin_index       (e_spin);

//     assert(configuration.get_changed_spin_indices().size()            == 0);
//     assert(configuration.get_changed_spin_indices_e_spin(e_UP).size() == 0);
//     assert(configuration.get_changed_spin_indices_e_spin(e_DN).size() == 0);
//     assert(configuration.assert_block_form(e_spin));
//     assert(first_shuffled_vertex_index >= first_non_interacting_vertex_index);

//     {
//       N.resize(configuration_size);

//       assert(   configuration_size == N.get_current_size());
//       assert(G0.get_current_size() == N.get_current_size());
//     }

//     {
//       //profiler_t profiler(concurrency, "TOOLS update_N_matrix b reset", __FILE__, __LINE__);
//       for(int j=0; j<configuration_size; j++){
//      if(configuration_e_spin[j].get_HS_spin() == HS_ZERO){
//        memset(&N(0,j), 0, sizeof(double)*configuration_size);


//        N(j,j) = double(1);
//      }
//       }
//     }

//     if(first_shuffled_vertex_index == configuration_size)
//       return;

//     if(first_non_interacting_vertex_index == 0)
//       return;

//     {// \sum_{l=0}^{l<vertex_index} exp_V_minus_one[l] * N(l,j) ---> exp_V_minus_one_times_N
//       //profiler_t profiler(concurrency, "TOOLS update_N_matrix c resize G0_times_exp_V_minus_one", __FILE__, __LINE__);
//       G0_times_exp_V_minus_one.resize_no_copy(std::pair<int,int>(configuration_size-first_shuffled_vertex_index, first_non_interacting_vertex_index));
//     }

//     {
//       for(int j=0; j<first_non_interacting_vertex_index; j++){

//      HS_spin_states_type HS_spin  = configuration_e_spin[j].get_HS_spin();
//      HS_field_sign_type  HS_field = configuration_e_spin[j].get_HS_field();

//      int spin_orbital         = configuration_e_spin[j].get_spin_orbital();
//      int spin_orbital_paired  = configuration_e_spin[j].get_paired_spin_orbital();

//      double  exp_V_minus_one_val = CV_obj.exp_V(spin_orbital, spin_orbital_paired, HS_spin, HS_field)-double(1);

//      double* G0_times_exp_V_minus_one_ptr = &G0_times_exp_V_minus_one(0,j);
//      memcpy(G0_times_exp_V_minus_one_ptr, &G0(first_shuffled_vertex_index,j), sizeof(double)*(configuration_size-first_shuffled_vertex_index));

//      for(int i=0; i<(configuration_size-first_shuffled_vertex_index); i++)
//        G0_times_exp_V_minus_one(i,j) *= exp_V_minus_one_val;
//       }
//     }

//     { // G0_exp_V_minus_one * N

//       //gemm_plan<double, BLAS_LIBRARY> dgemm_pln;
//       gemm_plan<double> dgemm_pln;
//       {
//      dgemm_pln.A   = &G0_times_exp_V_minus_one(0,0);
//      dgemm_pln.LDA = G0_times_exp_V_minus_one.get_global_size().first;

//      dgemm_pln.B   = &N(0,0);
//      dgemm_pln.LDB = N.get_global_size();

//      dgemm_pln.C   = &N(first_shuffled_vertex_index,0);
//      dgemm_pln.LDC = N.get_global_size();

//      dgemm_pln.M   = configuration_size - first_shuffled_vertex_index;
//      dgemm_pln.K   = first_non_interacting_vertex_index;
//      dgemm_pln.N   = first_non_interacting_vertex_index;

//      //dgemm_pln.execute_plan();
//      dgemm_pln.execute_plan(my_blas_library);
//       }
//     }
//   }


//   template<LIN_ALG::device_type device_t, typename parameters_type, class base_cluster_type>
//   template<class configuration_type, class vertex_vertex_matrix_type, class vertex_new_vertex_matrix_type>
//   void  N_TOOLS<device_t, parameters_type, base_cluster_type>::rebuild_N_matrix_via_Gamma_LU(configuration_type&            full_configuration,
//                                                                                           vertex_vertex_matrix_type&     N,
//                                                                                           vertex_vertex_matrix_type&     Gamma,
//                                                                                           vertex_new_vertex_matrix_type& G_precomputed,
//                                                                                           e_spin_states_type             e_spin)
//   {
//     profiler_t profiler(concurrency, "TOOLS rebuild_N_matrix_via_Gamma_LU", __FILE__, __LINE__);

//     int Gamma_size = Gamma.get_current_size();

//     if(Gamma_size == 0 )
//       return;

//     std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
//     int                                 configuration_size   = configuration_e_spin.size();

//     std::vector<HS_spin_states_type>& spin_values = full_configuration.get_changed_spin_values_e_spin(e_spin);
//     std::vector<int>&                 permutation = full_configuration.get_changed_spin_indices_e_spin(e_spin);

//     assert(spin_values.size()   ==     permutation.size());
//     assert(Gamma_size           == int(permutation.size()));
//     assert(N.get_current_size() == int(configuration_size));
//     assert(assert_that_there_are_no_Bennett_spins(full_configuration));

//     { // get the rows of N corresponding to the new spins => N_new_spins
//       profiler_t profiler(concurrency, "TOOLS rebuild_N_matrix_via_Gamma_LU a resize N && copy rows", __FILE__, __LINE__);

//       N_new_spins.resize_no_copy(std::pair<int,int>(Gamma_size, configuration_size));

//       copy_rows_N(permutation, N);
// //       for(int j=0;j<configuration_size;j++)
// //   for(int i=0;i<Gamma_size;i++)
// //     N_new_spins(i,j) = N(permutation[i],j);
//     }

//     { // get the columns of G corresponding to the new spins => G_new_spins
//       profiler_t profiler(concurrency, "TOOLS rebuild_N_matrix_via_Gamma_LU b resize G && copy cols", __FILE__, __LINE__);

//       G.resize_no_copy(std::pair<int,int>(configuration_size, Gamma_size));

//       compute_G_changed_vertices_to_all_vertex(N, G_precomputed, full_configuration, e_spin);
//     }

//     { // Gamma_LU * X = N(p_k,:) --> X = Gamma_inv_times_N_new_spins ==> (stored in N_new_spins)
//       profiler_t profiler(concurrency, "TOOLS rebuild_N_matrix_via_Gamma_LU c LU-solve", __FILE__, __LINE__);

//       getrs_plan<double> getrs_pln(Gamma_size, configuration_size, Gamma.get_global_size(), N_new_spins.get_global_size().first);
//       getrs_pln.Matrix_A = &(Gamma(0,0));
//       getrs_pln.Matrix_B = &(N_new_spins(0,0));
//       getrs_pln.execute_plan(my_lapack_library);
//     }

//     // do N - G*Gamma_inv_times_N_new_spins --> N  || DGEMM --> work-horse
//     {
//       profiler_t profiler(concurrency, "TOOLS rebuild_N_matrix_via_Gamma_LU d dgemm", __FILE__, __LINE__);

//       //gemm_plan<double, BLAS_LIBRARY>   gemm_pln(configuration_size);
//       gemm_plan<double>   gemm_pln(configuration_size);
//       {
//      gemm_pln.M = configuration_size;
//      gemm_pln.N = configuration_size;
//      gemm_pln.K = Gamma_size;

//      gemm_pln.alpha = -1.;
//      gemm_pln.beta  =  1.;

//      gemm_pln.A = &(G(0,0));
//      gemm_pln.B = &(N_new_spins(0,0));
//      gemm_pln.C = &(N(0,0));

//      gemm_pln.LDA = G          .get_global_size().first;
//      gemm_pln.LDB = N_new_spins.get_global_size().first;
//      gemm_pln.LDC = N          .get_global_size();

//      //gemm_pln.execute_plan();
//      gemm_pln.execute_plan(my_blas_library);
//       }
//     }

//     {
//       profiler_t profiler(concurrency, "TOOLS rebuild_N_matrix_via_Gamma_LU e rescale", __FILE__, __LINE__);

//       // do N*D_i --> N ( = final N !)
//       std::vector<double> d_inv(permutation.size());

//       for(int i=0;i<int(permutation.size());i++)
//      {
//        HS_spin_states_type old_HS_spin   = configuration_e_spin[permutation[i]].get_HS_spin();
//        HS_spin_states_type new_HS_spin   = spin_values[i];

//        HS_field_sign       HS_field_sign = configuration_e_spin[permutation[i]].get_HS_field();

//        int spin_orbital                  = configuration_e_spin[permutation[i]].get_spin_orbital();
//        int spin_orbital_paired           = configuration_e_spin[permutation[i]].get_paired_spin_orbital();

//        if(old_HS_spin == HS_ZERO)
//          {
//            double exp_delta_V = CV_obj.exp_delta_V(spin_orbital, spin_orbital_paired, new_HS_spin, old_HS_spin, HS_field_sign);

//            double phani_gamma = exp_delta_V-1.;

//            d_inv[i] = 1./(1. + phani_gamma);
//          }
//        else
//          {
//            d_inv[i] = 1./N(permutation[i],permutation[i]);
//          }
//      }

//       scale_rows_N(permutation, d_inv, N);
// //       for(int j=0;j<N.get_current_size();++j){
// //   for(int i=0;i< int(permutation.size());++i){
// //     N(permutation[i],j) *= d_inv[i];
// //   }
// //       }
//     }
//   }

//   template<LIN_ALG::device_type device_t, typename parameters_type, class base_cluster_type>
//   template<class vertex_vertex_matrix_type>
//   inline void N_TOOLS<device_t, parameters_type, base_cluster_type>::copy_rows_N(std::vector<int>&          permutation,
//                                                                               vertex_vertex_matrix_type& N)
//   {
// //     int Gamma_size         = permutation.size();
// //     int configuration_size = N.get_current_size();

// //     for(int j=0;j<configuration_size;j++)
// //       for(int i=0;i<Gamma_size;i++)
// //   N_new_spins(i,j) = N(permutation[i],j);

//     int length = N.get_current_size();

//     double* N_ptr           = &N(0, 0);
//     double* N_new_spins_ptr = &N_new_spins(0, 0);

//     int N_LD     = N.get_global_size();
//     int N_new_LD = N_new_spins.get_global_size().first;

//     for(size_t i=0;i<permutation.size();i++)
//       copy_plan::execute(length, &N_ptr[permutation[i]], N_LD, &N_new_spins_ptr[i], N_new_LD);
//   }

//   template<LIN_ALG::device_type device_t, typename parameters_type, class base_cluster_type>
//   template<class vertex_vertex_matrix_type>
//   inline void N_TOOLS<device_t, parameters_type, base_cluster_type>::scale_rows_N(std::vector<int>&              permutation,
//                                                                                std::vector<double>&       d_inv,
//                                                                                vertex_vertex_matrix_type& N)
//   {
// //     for(int j=0;j<N.get_current_size();++j)
// //       for(int i=0;i< int(permutation.size());++i)
// //   N(permutation[i],j) *= d_inv[i];

//     int length = N.get_current_size();
//     int N_LD   = N.get_global_size();

//     double* N_ptr = &N(0, 0);

//     for(size_t i=0;i<permutation.size();i++)
//       scale_plan::execute(length, d_inv[i], &N_ptr[permutation[i]], N_LD);
//   }

//   template<LIN_ALG::device_type device_t, typename parameters_type, class base_cluster_type>
//   template<class configuration_type, class vertex_vertex_matrix_type, class vertex_new_vertex_matrix_type>
//   void N_TOOLS<device_t, parameters_type, base_cluster_type>::compute_G_changed_vertices_to_all_vertex(vertex_vertex_matrix_type&     N,
//                                                                                                     vertex_new_vertex_matrix_type& G_precomputed,
//                                                                                                     configuration_type&            full_configuration,
//                                                                                                     e_spin_states_type             e_spin)
//   {
//     std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

//     std::vector<int>&                   permutation = full_configuration.get_changed_spin_indices_e_spin(e_spin);

//     int                                 configuration_size   = configuration_e_spin.size();
//     int                                 Gamma_size           = permutation.size();

//     int vertex_index=0;
//     while(vertex_index < configuration_size && configuration_e_spin[vertex_index].get_HS_spin() != HS_ZERO)
//       vertex_index++;

//     for(int j=0; j<Gamma_size; j++)
//       {
//      int configuration_index_e_spin = permutation[j];

//      HS_spin_states_type HS_spin_2 = configuration_e_spin[configuration_index_e_spin].get_HS_spin();

//      if(HS_spin_2 != HS_ZERO)
//        {
//          // Gij = (Nij e^Vj - delta_ij)/(e^Vj -1) (eqn 33)

//          HS_spin_states_type HS_spin  = configuration_e_spin[configuration_index_e_spin].get_HS_spin();
//          HS_field_sign_type  HS_field = configuration_e_spin[configuration_index_e_spin].get_HS_field();

//          int spin_orbital         = configuration_e_spin[configuration_index_e_spin].get_spin_orbital();
//          int spin_orbital_paired  = configuration_e_spin[configuration_index_e_spin].get_paired_spin_orbital();

//          double exp_V = CV_obj./*CV::*/exp_V(spin_orbital, spin_orbital_paired, HS_spin, HS_field);

//          //double numerator   = exp_V - delta;
//          double denumerator = exp_V-1.;
//          double delta;

//          for(int i=0; i<configuration_size; i++){
//            delta = i==configuration_index_e_spin? 1.:0;

//            G(i,j) = N(i,configuration_index_e_spin) * (exp_V - delta)/denumerator;
//          }
//        }

//      if(HS_spin_2 == HS_ZERO)
//        {
//          // G_i_j = \sum_l N(configuration_e_spin_index_i,l) * G0(l,configuration_e_spin_index_j);
//          memcpy(&G(0,j), &G_precomputed(0,configuration_index_e_spin-vertex_index), sizeof(double)*configuration_size);
//        }
//       }

//   }

/*
  template<LIN_ALG::device_type device_t>
  class N_TOOLS_data
  {};

  template<>
  class N_TOOLS_data<LIN_ALG::CPU>
  {
  public :

  N_TOOLS_data() {}
  ~N_TOOLS_data() {}

  double* get_device_ptr(LIN_ALG::vector<double, LIN_ALG::CPU>& v){
  return v.get_ptr();
  }
  };

  template<>
  class N_TOOLS_data<LIN_ALG::GPU>
  {
  public :

  N_TOOLS_data():
  id(128),
  permutation(128),
  tmp(128)
  {
  std::vector<int> id_tmp(128);
  for(int l=0; l<128; ++l)
  id_tmp[l] = l;

  id.set(id_tmp);
  }

  ~N_TOOLS_data()
  {}

  void set_permutation(std::vector<int>& p){
  permutation.set(p);
  }

  int* get_permutation(){
  return permutation.get_ptr();
  }

  double* get_device_ptr(LIN_ALG::vector<double, LIN_ALG::CPU>& v){

  tmp.resize(v.get_current_size());

  LIN_ALG::COPY_FROM<LIN_ALG::CPU, LIN_ALG::GPU>::execute(v.get_ptr(), tmp.get_ptr(), v.get_current_size());

  return tmp.get_ptr();
  }

  void copy_rows(LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
  LIN_ALG::matrix<double, LIN_ALG::GPU>& N_new_spins);

  void copy_cols(LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
  LIN_ALG::matrix<double, LIN_ALG::GPU>& G,
  LIN_ALG::matrix<double, LIN_ALG::GPU>& G_new_spins);

  void scale_rows(std::vector<double>&                   alpha,
  LIN_ALG::matrix<double, LIN_ALG::GPU>& G_new_spins);

  void scale_cols(std::vector<double>&                   alpha,
  LIN_ALG::matrix<double, LIN_ALG::GPU>& G_new_spins);

  private:

  LIN_ALG::vector<int, LIN_ALG::GPU> id;
  LIN_ALG::vector<int, LIN_ALG::GPU> permutation;

  LIN_ALG::vector<double, LIN_ALG::GPU> tmp;
  LIN_ALG::vector<double, LIN_ALG::GPU> alpha;
  };
*/
