//-*-C++-*-

#ifndef DCA_QMCI_ACCUMULATOR_NONLOCAL_G_H
#define DCA_QMCI_ACCUMULATOR_NONLOCAL_G_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  namespace QMCI
  {
    namespace CT_AUX_ACCUMULATION
    {

      /*!
       *  \class   accumulator_nonlocal_G
       *  \ingroup CT-AUX-ACCUMULATOR
       *
       *  \brief   This class computes the nonlocal single-particle Greens-function \f$G^{I}(k_1,k_2)\f$
       *  \date    12-07-2011
       *  \author  Peter Staar
       *  \version 1.0
       */
      template<class parameters_type, class MOMS_type>
      class accumulator_nonlocal_G
      {

        typedef r_DCA r_dmn_t;
        typedef k_DCA k_dmn_t;

        typedef vertex_singleton        vertex_singleton_type;

        typedef typename parameters_type::profiler_type    profiler_t;
        typedef typename parameters_type::concurrency_type concurrency_type;

        typedef typename parameters_type::MC_measurement_scalar_type scalar_type;

        typedef typename parameters_type::G4_w1_dmn_t w1_dmn_t;
        typedef typename parameters_type::G4_w2_dmn_t w2_dmn_t;

        typedef dmn_6<b,b,r_dmn_t,r_dmn_t,w1_dmn_t,w2_dmn_t> b_b_r_r_w_w_dmn_t;
        typedef dmn_6<b,b,k_dmn_t,r_dmn_t,w1_dmn_t,w2_dmn_t> b_b_k_r_w_w_dmn_t;
        typedef dmn_6<b,b,k_dmn_t,k_dmn_t,w1_dmn_t,w2_dmn_t> b_b_k_k_w_w_dmn_t;

        typedef LIN_ALG::matrix<double, LIN_ALG::CPU> vertex_vertex_matrix_type;

      public:

        accumulator_nonlocal_G(parameters_type& parameters_ref,
                               MOMS_type&       MOMS_ref,
                               int              id);

        ~accumulator_nonlocal_G();

        double get_GFLOP();

        void initialize();

        void finalize();

        /*!
         * \brief Compute the nonlocal single particle Greens-function \f$G^{I}_{\sigma}(k_1, k_2)\f$, from the M-matrices.
         */
        void execute(std::vector<vertex_singleton_type>& HS_configuration_e_UP, vertex_vertex_matrix_type& M_e_UP,
                     std::vector<vertex_singleton_type>& HS_configuration_e_DN, vertex_vertex_matrix_type& M_e_DN);

        /*!
         * \brief Get a reference to \f$G^{I}_{\uparrow}(k_1, k_2)\f$.
         */
        FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& get_G_k_k_w_w_e_UP()
        { return G2_k_k_w_w_e_UP; }

        /*!
         * \brief Get a reference to \f$G^{I}_{\downarrow}(k_1, k_2)\f$.
         */
        FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& get_G_k_k_w_w_e_DN()
        { return G2_k_k_w_w_e_DN; }

        /*!
         * \brief print help-functions.
         */
        template<class stream_type>
        void to_JSON(stream_type& ss);

      private:

        void FT_M_v_v_2_M_k_k_w_w(std::vector<vertex_singleton_type>& HS_configuration_e_spin,
                                  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w,
                                  vertex_vertex_matrix_type& M_e_spin);

        void FT_M_v_v_2_M_k_k_w_w_test();

        void compute_G2_k_k_w_w(FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w,
                                FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w,
                                e_spin_states_type e_spin);

        void compute_G2_k_k_w_w_atomic(std::complex<double>* G0_k_w_l,
                                       std::complex<scalar_type>* M_k_w,
                                       std::complex<double>* G0_k_w_r,
                                       std::complex<scalar_type>* G2_k_k_w_w);

        void compute_G2_k_k_w_w_atomic_1(std::complex<double>* G0_k_w_l,
                                         std::complex<scalar_type>* M_k_w,
                                         std::complex<double>* G0_k_w_r,
                                         std::complex<scalar_type>* G2_k_k_w_w);

        void compute_G2_k_k_w_w_atomic_default(std::complex<double>* G0_k_w_l,
                                               std::complex<scalar_type>* M_k_w,
                                               std::complex<double>* G0_k_w_r,
                                               std::complex<scalar_type>* G2_k_k_w_w);

      private:

        parameters_type&          parameters;
        MOMS_type&                MOMS;
        concurrency_type&         concurrency;

        int                       thread_id;

        double GFLOP;

        scalar_type beta;

        int b_dmn_size;
        int H_dmn_size;

        std::vector<int> corresponding_w1_index;
        std::vector<int> corresponding_w2_index;
	cached_nft<2, scalar_type, r_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED_POS> nft_obj;

        FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t> G2_k_k_w_w_e_UP;
        FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t> G2_k_k_w_w_e_DN;

        FUNC_LIB::function<std::complex<scalar_type>, b_b_r_r_w_w_dmn_t> M_r_r_w_w;
        FUNC_LIB::function<std::complex<scalar_type>, b_b_k_r_w_w_dmn_t> M_k_r_w_w;
        FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t> M_k_k_w_w;
      };

      template<class parameters_type, class MOMS_type>
      accumulator_nonlocal_G<parameters_type, MOMS_type>::accumulator_nonlocal_G(parameters_type&          parameters_ref,
                                                                                 MOMS_type&                MOMS_ref,
                                                                                 int                       id):
        parameters(parameters_ref),
        MOMS(MOMS_ref),
        concurrency(parameters.get_concurrency()),

        thread_id(id),

        GFLOP(0),

        beta(parameters.get_beta()),

        b_dmn_size(b::dmn_size()),
        H_dmn_size(b::dmn_size()*s::dmn_size()),

        corresponding_w1_index(w1_dmn_t::parameter_type::get_corresponding_frequency_domain_index()),
        corresponding_w2_index(w2_dmn_t::parameter_type::get_corresponding_frequency_domain_index()),

        nft_obj(w_VERTEX_EXTENDED::parameter_type::get_elements()),

        G2_k_k_w_w_e_UP("G2_k_k_w_w_e_UP"),
        G2_k_k_w_w_e_DN("G2_k_k_w_w_e_DN"),

        M_r_r_w_w("M_r_r_w_w_tmp"),
        M_k_r_w_w("M_k_r_w_w_tmp"),
        M_k_k_w_w("M_k_k_w_w_tmp")
      {}

      template<class parameters_type, class MOMS_type>
      accumulator_nonlocal_G<parameters_type, MOMS_type>::~accumulator_nonlocal_G()
      {}

      template<class parameters_type, class MOMS_type>
      double accumulator_nonlocal_G<parameters_type, MOMS_type>::get_GFLOP()
      {
        double result = GFLOP;
        GFLOP        = 0;

        return result;
      }

      template<class parameters_type, class MOMS_type>
      void accumulator_nonlocal_G<parameters_type, MOMS_type>::initialize()
      {
        beta = parameters.get_beta();

        b_dmn_size = b::dmn_size();
        H_dmn_size = b::dmn_size()*s::dmn_size();
      }

      template<class parameters_type, class MOMS_type>
      void accumulator_nonlocal_G<parameters_type, MOMS_type>::finalize()
      {}

      template<class parameters_type, class MOMS_type>
      void accumulator_nonlocal_G<parameters_type, MOMS_type>::execute(std::vector<vertex_singleton_type>& HS_configuration_e_UP, vertex_vertex_matrix_type& M_e_UP,
                                                                       std::vector<vertex_singleton_type>& HS_configuration_e_DN, vertex_vertex_matrix_type& M_e_DN)
      {
        profiler_t profiler("compute nonlocal-G from M-matrix", "CT-AUX accumulator", __LINE__, thread_id);

        //e_UP          
	FT_M_v_v_2_M_k_k_w_w(HS_configuration_e_UP, M_k_k_w_w, M_e_UP);
	
	//profiler_t profiler("(B)\tb. c. --> nonlocal-G-compute       ", __FILE__, __LINE__);
	compute_G2_k_k_w_w(M_k_k_w_w, G2_k_k_w_w_e_UP, e_UP);
        
	//e_DN
	//profiler_t profiler("(B)\tb. c. --> nonlocal-G-FT-M          ", __FILE__, __LINE__);
	FT_M_v_v_2_M_k_k_w_w(HS_configuration_e_DN, M_k_k_w_w, M_e_DN);
          
	//profiler_t profiler("(B)\tb. c. --> nonlocal-G-compute       ", __FILE__, __LINE__);
	compute_G2_k_k_w_w(M_k_k_w_w, G2_k_k_w_w_e_DN, e_DN);
      }

      template<class parameters_type, class MOMS_type>
      void accumulator_nonlocal_G<parameters_type, MOMS_type>::FT_M_v_v_2_M_k_k_w_w(std::vector<vertex_singleton_type>&                     configuration_e_spin,
                                                                                    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w_e_spin,
                                                                                    vertex_vertex_matrix_type&                              M_e_spin)
      {
        {
          profiler_t profiler("nonlocal-G-FT t --> w", "CT-AUX accumulator", __LINE__, thread_id);
          GFLOP += nft_obj.execute(configuration_e_spin, M_e_spin, M_r_r_w_w);
        }

        {
          profiler_t profiler("nonlocal-G-FT r --> k", "CT-AUX accumulator", __LINE__, thread_id);

          math_algorithms::functional_transforms::TRANSFORM<r_dmn_t, k_dmn_t>::execute(M_r_r_w_w, M_k_r_w_w);
          math_algorithms::functional_transforms::TRANSFORM<r_dmn_t, k_dmn_t>::execute(M_k_r_w_w, M_k_k_w_w_e_spin);

          scalar_type  Nc        = scalar_type(k_DCA::dmn_size());
          scalar_type one_div_Nc = 1./Nc;
          M_k_k_w_w_e_spin *= one_div_Nc;

          GFLOP += M_r_r_w_w.size()*double(Nc)*(1.e-9);
        }
      }

    
      template<class parameters_type, class MOMS_type>
      void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w(FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w,
                                                                                  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w,
                                                                                  e_spin_states_type e_spin)
      {
        std::complex<double>      *G0_k_w_l_ptr , *G0_k_w_r_ptr;
        std::complex<scalar_type> *M_k_k_w_w_ptr, *G_k_k_w_w_ptr;

	int spin_index = electron_spin_domain::to_coordinate(e_spin);

        for(int i=0; i<G2_k_k_w_w.size(); i++)
          G2_k_k_w_w(i) = 0;

        for(int w2=0; w2<w2_dmn_t::dmn_size(); w2++){

          int corr_w2 = corresponding_w2_index[w2];

          for(int w1=0; w1<w1_dmn_t::dmn_size(); w1++){

            int corr_w1 = corresponding_w1_index[w1];

            for(int k2=0; k2<k_dmn_t::dmn_size(); k2++){
              for(int k1=0; k1<k_dmn_t::dmn_size(); k1++){

                G0_k_w_l_ptr  = &MOMS.G0_k_w_cluster_excluded(0, spin_index, 0, spin_index, k1, corr_w1);
                M_k_k_w_w_ptr = &M_k_k_w_w(0, 0, k1, k2, w1, w2);
                G0_k_w_r_ptr  = &MOMS.G0_k_w_cluster_excluded(0, spin_index, 0, spin_index, k2, corr_w2);

                G_k_k_w_w_ptr = &G2_k_k_w_w(0, 0, k1, k2, w1, w2);

                compute_G2_k_k_w_w_atomic(G0_k_w_l_ptr, M_k_k_w_w_ptr, G0_k_w_r_ptr, G_k_k_w_w_ptr);
              }
            }
          }
        }
      }

      template<class parameters_type, class MOMS_type>
      inline void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w_atomic(std::complex<double>* G0_k_w_l_ptr,
                                                                                                std::complex<scalar_type>* M_k_k_w_w_ptr,
                                                                                                std::complex<double>* G0_k_w_r_ptr,
                                                                                                std::complex<scalar_type>* G2_k_k_w_w_ptr)
      {
        switch(b_dmn_size)
          {
          case 1:
            compute_G2_k_k_w_w_atomic_1(G0_k_w_l_ptr, M_k_k_w_w_ptr, G0_k_w_r_ptr, G2_k_k_w_w_ptr);
            break;

          default:
            compute_G2_k_k_w_w_atomic_default(G0_k_w_l_ptr, M_k_k_w_w_ptr, G0_k_w_r_ptr, G2_k_k_w_w_ptr);
          }
      }

      template<class parameters_type, class MOMS_type>
      inline void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w_atomic_1(std::complex<double>* G0_k_w_l_ptr,
                                                                                                  std::complex<scalar_type>* M_k_k_w_w_ptr,
                                                                                                  std::complex<double>* G0_k_w_r_ptr,
                                                                                                  std::complex<scalar_type>* G2_k_k_w_w_ptr)
      {
        G2_k_k_w_w_ptr[0] -= (std::complex<scalar_type>(G0_k_w_l_ptr[0]) * M_k_k_w_w_ptr[0] * std::complex<scalar_type>(G0_k_w_r_ptr[0]));

        if(G0_k_w_l_ptr == G0_k_w_r_ptr)
          G2_k_k_w_w_ptr[0] += std::complex<scalar_type>(G0_k_w_l_ptr[0])*beta;
      }

      template<class parameters_type, class MOMS_type>
      inline void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w_atomic_default(std::complex<double>* G0_k_w_l_ptr,
                                                                                                        std::complex<scalar_type>* M_k_k_w_w_ptr,
                                                                                                        std::complex<double>* G0_k_w_r_ptr,
                                                                                                        std::complex<scalar_type>* G2_k_k_w_w_ptr)
      {
        for(int b1=0; b1<b_dmn_size; b1++){
          for(int b2=0; b2<b_dmn_size; b2++){

            for(int l1=0; l1<b_dmn_size; l1++)
              for(int l2=0; l2<b_dmn_size; l2++)
                G2_k_k_w_w_ptr[b1+b2*b_dmn_size] -= (std::complex<scalar_type>(G0_k_w_l_ptr[b1+l1*H_dmn_size])
                                                     * M_k_k_w_w_ptr[l1+l2*b_dmn_size] *
                                                     std::complex<scalar_type>(G0_k_w_r_ptr[l2+b2*H_dmn_size]));

            if(G0_k_w_l_ptr == G0_k_w_r_ptr)
              G2_k_k_w_w_ptr[b1+b2*b_dmn_size] += std::complex<scalar_type>(G0_k_w_l_ptr[b1+b2*H_dmn_size])*beta;
          }
        }
      }

    }

  }

}

#endif
